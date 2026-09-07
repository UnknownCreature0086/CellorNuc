#!/usr/bin/env python3
"""Run CellorNucEM on quiz dataset A and summarize Q1a.

The fitted mixture is unsupervised with respect to ``suspension_type``.
That metadata is passed to CellorNucEM only so the final labels distinguish
agreement ("Typical") from cross-modality calls. We additionally derive a
plain cell-like/nucleus-like prediction from ``p_cell`` so accuracy and
coverage can be evaluated without relying on the metadata-aware label names.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import average_precision_score, roc_auc_score


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scanpyplus import Scanpyplus  # noqa: E402


DEFAULT_INPUT = PROJECT_ROOT / "data/mixed_mBDRC_5pCells_multiomeNuclei_6k.h5ad"
DEFAULT_RESULTS = PROJECT_ROOT / "results/q1"
DEFAULT_FIGURE = PROJECT_ROOT / "figures/q1_cellornucem_diagnostics.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--min-counts", type=int, default=10)
    parser.add_argument("--gamma-hi", type=float, default=0.95)
    parser.add_argument("--gamma-lo", type=float, default=0.05)
    parser.add_argument(
        "--skip-h5ad",
        action="store_true",
        help="Do not save the annotated AnnData result.",
    )
    return parser.parse_args()


def component_mean(ab: tuple[float, float] | list[float]) -> float:
    """Return the mean of a beta distribution parameterized by (a, b)."""
    a, b = ab
    return float(a / (a + b))


def records_by_index(table: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Convert a crosstab to a JSON-friendly nested mapping."""
    return {
        str(index): {str(column): int(value) for column, value in row.items()}
        for index, row in table.iterrows()
    }


def make_fit_h5ad_serializable(value: object) -> object:
    """Convert tuple-valued EM parameters into AnnData-writable arrays."""
    if isinstance(value, dict):
        return {key: make_fit_h5ad_serializable(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return np.asarray(value)
    return value


def plot_diagnostics(
    adata: sc.AnnData,
    prediction: pd.Categorical,
    output_path: Path,
    min_counts: int,
    gamma_lo: float,
    gamma_hi: float,
) -> None:
    """Plot known assay, independent calls, and signature separation."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth = adata.obs["suspension_type"].astype(str).to_numpy()
    pred = np.asarray(prediction).astype(str)
    umap = adata.obsm["X_umap"]

    truth_colors = {"cell": "#e76f51", "nucleus": "#457b9d"}
    pred_colors = {
        "Cell-like": "#e76f51",
        "Nucleus-like": "#457b9d",
        "Neutral": "#b7b7b7",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    for label, color in truth_colors.items():
        mask = truth == label
        axes[0, 0].scatter(
            umap[mask, 0], umap[mask, 1], s=3, alpha=0.65,
            linewidths=0, color=color, label=label,
        )
    axes[0, 0].set_title("Known suspension type")
    axes[0, 0].legend(markerscale=4, frameon=False)

    for label, color in pred_colors.items():
        mask = pred == label
        axes[0, 1].scatter(
            umap[mask, 0], umap[mask, 1], s=3, alpha=0.65,
            linewidths=0, color=color, label=label,
        )
    axes[0, 1].set_title("CellorNucEM prediction (truth not used in fit)")
    axes[0, 1].legend(markerscale=4, frameon=False)

    informative = adata.obs["modality_counts"].to_numpy() >= min_counts
    for label, color in truth_colors.items():
        values = adata.obs.loc[informative & (truth == label), "sc_frac"]
        axes[1, 0].hist(
            values, bins=np.linspace(0, 1, 51), density=True,
            histtype="step", linewidth=2, color=color, label=label,
        )
    axes[1, 0].set(
        title="Observed cell-signature fraction (informative droplets)",
        xlabel="sc-signature counts / all signature counts",
        ylabel="Density",
        xlim=(0, 1),
    )
    axes[1, 0].legend(frameon=False)

    for label, color in truth_colors.items():
        values = adata.obs.loc[truth == label, "p_cell"]
        axes[1, 1].hist(
            values, bins=np.linspace(0, 1, 51), density=True,
            histtype="step", linewidth=2, color=color, label=label,
        )
    axes[1, 1].axvline(gamma_lo, color="black", linestyle="--", linewidth=1)
    axes[1, 1].axvline(gamma_hi, color="black", linestyle="--", linewidth=1)
    axes[1, 1].set(
        title="Posterior probability of the cell-like component",
        xlabel="p_cell",
        ylabel="Density",
        xlim=(0, 1),
    )
    axes[1, 1].legend(frameon=False)

    for ax in axes[0]:
        ax.set(xlabel="UMAP1", ylabel="UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("Q1a: CellorNucEM on mixed dataset A", fontsize=15)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    if "suspension_type" not in adata.obs.columns:
        raise ValueError("Missing required obs column: suspension_type")
    if "counts" not in adata.layers:
        raise ValueError("Dataset has no raw-count layer named 'counts'.")

    sc_genes, sn_genes = Scanpyplus.LoadGeneSignatures()
    sc_found = [gene for gene in sc_genes if gene in adata.var_names]
    sn_found = [gene for gene in sn_genes if gene in adata.var_names]

    # This follows the professor-provided mixed-data demo: a global fit on raw
    # counts with conservative 0.95/0.05 posterior thresholds. The exact
    # metadata values in this file are "cell" and "nucleus", not "sc"/"sn".
    result = Scanpyplus.CellorNucEM(
        adata,
        sc_genes,
        sn_genes,
        layer="counts",
        assay_col="suspension_type",
        sc_label="cell",
        sn_label="nucleus",
        min_counts=args.min_counts,
        gamma_hi=args.gamma_hi,
        gamma_lo=args.gamma_lo,
        print_mode=False,
        copy=True,
    )

    truth = result.obs["suspension_type"].astype(str)
    truth_values = set(truth.unique())
    if truth_values != {"cell", "nucleus"}:
        raise ValueError(
            "Expected suspension_type values {'cell', 'nucleus'}, "
            f"found {sorted(truth_values)}"
        )
    classification = result.obs["modality_classification"].astype(str)
    enough = result.obs["modality_counts"].to_numpy() >= args.min_counts
    p_cell = result.obs["p_cell"].to_numpy()

    predicted = np.full(result.n_obs, "Neutral", dtype=object)
    predicted[enough & (p_cell > args.gamma_hi)] = "Cell-like"
    predicted[enough & (p_cell < args.gamma_lo)] = "Nucleus-like"
    prediction = pd.Categorical(
        predicted, categories=["Cell-like", "Nucleus-like", "Neutral"]
    )
    result.obs["predicted_modality"] = prediction

    class_order = [
        "Typical Cell (SC)",
        "Nucleus-like Cell (SC)",
        "Typical Nucleus (SN)",
        "Cell-like Nucleus (SN)",
        "Neutral",
    ]
    classification_counts = classification.value_counts().reindex(
        class_order, fill_value=0
    )
    assay_by_classification = pd.crosstab(truth, classification).reindex(
        index=["cell", "nucleus"], columns=class_order, fill_value=0
    )
    assay_by_prediction = pd.crosstab(truth, prediction).reindex(
        index=["cell", "nucleus"],
        columns=["Cell-like", "Nucleus-like", "Neutral"],
        fill_value=0,
    )

    confident = predicted != "Neutral"
    correct = ((truth.to_numpy() == "cell") & (predicted == "Cell-like")) | (
        (truth.to_numpy() == "nucleus") & (predicted == "Nucleus-like")
    )
    n_confident = int(confident.sum())
    cell_mask = truth.to_numpy() == "cell"
    nucleus_mask = truth.to_numpy() == "nucleus"
    n_cell_confident = int((confident & cell_mask).sum())
    n_nucleus_confident = int((confident & nucleus_mask).sum())
    y_cell = cell_mask.astype(int)
    fit = result.uns["modality_em"]["global"]
    cell_mean = component_mean(fit["ab_cell"])
    nucleus_mean = component_mean(fit["ab_nuc"])

    cell_total = int((truth == "cell").sum())
    nucleus_total = int((truth == "nucleus").sum())
    cell_as_nucleus = int(
        ((truth == "cell") & (classification == "Nucleus-like Cell (SC)")).sum()
    )
    nucleus_as_cell = int(
        ((truth == "nucleus") & (classification == "Cell-like Nucleus (SN)")).sum()
    )

    summary = {
        "input": str(args.input.resolve()),
        "shape": {"n_obs": int(result.n_obs), "n_vars": int(result.n_vars)},
        "parameters": {
            "fit": "global",
            "layer": "counts",
            "assay_col": "suspension_type",
            "sc_label": "cell",
            "sn_label": "nucleus",
            "min_counts": args.min_counts,
            "gamma_hi": args.gamma_hi,
            "gamma_lo": args.gamma_lo,
        },
        "signature_overlap": {
            "sc_found": sc_found,
            "sc_missing": sorted(set(sc_genes).difference(sc_found)),
            "sn_found": sn_found,
            "sn_missing": sorted(set(sn_genes).difference(sn_found)),
        },
        "fit_diagnostics": {
            "cell_component_mean_sc_fraction": cell_mean,
            "nucleus_component_mean_sc_fraction": nucleus_mean,
            "component_separation": cell_mean - nucleus_mean,
            "lambda_cell": float(fit["lam_cell"]),
            "bic_1_component": float(fit["bic_1"]),
            "bic_2_component": float(fit["bic_2"]),
            "delta_bic_1_minus_2": float(fit["bic_1"] - fit["bic_2"]),
        },
        "evidence": {
            "informative_droplets": int(enough.sum()),
            "low_signature_count_droplets": int((~enough).sum()),
            "informative_but_intermediate_posterior": int(
                ((predicted == "Neutral") & enough).sum()
            ),
            "low_signature_counts_by_assay": {
                "cell": int(((~enough) & cell_mask).sum()),
                "nucleus": int(((~enough) & nucleus_mask).sum()),
            },
        },
        "classification_counts": {
            str(key): int(value) for key, value in classification_counts.items()
        },
        "assay_by_classification": records_by_index(assay_by_classification),
        "assay_by_prediction": records_by_index(assay_by_prediction),
        "performance": {
            "confident_calls": n_confident,
            "coverage": n_confident / result.n_obs,
            "correct_confident_calls": int((correct & confident).sum()),
            "accuracy_among_confident_calls": (
                float((correct & confident).sum() / n_confident)
                if n_confident
                else None
            ),
            "correct_fraction_of_all_droplets": float(correct.sum() / result.n_obs),
            "accuracy_if_neutral_is_incorrect_among_informative": float(
                correct[enough].sum() / enough.sum()
            ),
            "coverage_among_informative": float(confident.sum() / enough.sum()),
            "auroc_on_informative_droplets": float(
                roc_auc_score(y_cell[enough], p_cell[enough])
            ),
            "average_precision_on_informative_droplets": float(
                average_precision_score(y_cell[enough], p_cell[enough])
            ),
            "balanced_accuracy_among_confident_calls": float(
                0.5
                * (
                    (n_cell_confident - cell_as_nucleus) / n_cell_confident
                    + (n_nucleus_confident - nucleus_as_cell) / n_nucleus_confident
                )
            ),
            "cell_called_nucleus_like": cell_as_nucleus,
            "cell_called_nucleus_like_rate": cell_as_nucleus / cell_total,
            "cell_confident_call_coverage": n_cell_confident / cell_total,
            "cell_accuracy_among_confident_calls": (
                (n_cell_confident - cell_as_nucleus) / n_cell_confident
            ),
            "nucleus_called_cell_like": nucleus_as_cell,
            "nucleus_called_cell_like_rate": nucleus_as_cell / nucleus_total,
            "nucleus_confident_call_coverage": n_nucleus_confident / nucleus_total,
            "nucleus_accuracy_among_confident_calls": (
                (n_nucleus_confident - nucleus_as_cell) / n_nucleus_confident
            ),
            # These literal combinations cannot be created by CellorNucEM when
            # assay_col is supplied; the two metrics above measure disagreement.
            "cell_labeled_literal_typical_nucleus": int(
                ((truth == "cell") & (classification == "Typical Nucleus (SN)")).sum()
            ),
            "nucleus_labeled_literal_typical_cell": int(
                ((truth == "nucleus") & (classification == "Typical Cell (SC)")).sum()
            ),
        },
    }

    classification_counts.rename("n_droplets").to_csv(
        args.results_dir / "classification_counts.csv"
    )
    assay_by_classification.to_csv(args.results_dir / "assay_by_classification.csv")
    assay_by_prediction.to_csv(args.results_dir / "assay_by_prediction.csv")
    result.obs[
        [
            "suspension_type",
            "cell_type",
            "modality_counts",
            "sc_frac",
            "p_cell",
            "predicted_modality",
            "modality_classification",
        ]
    ].to_csv(args.results_dir / "droplet_classifications.csv.gz", compression="gzip")
    with (args.results_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    plot_diagnostics(
        result,
        prediction,
        args.figure,
        args.min_counts,
        args.gamma_lo,
        args.gamma_hi,
    )
    if not args.skip_h5ad:
        # CellorNucEM stores (a, b) parameters as tuples, which AnnData's HDF5
        # writer cannot serialize. This changes only their container type.
        result.uns["modality_em"] = make_fit_h5ad_serializable(
            result.uns["modality_em"]
        )
        result.write_h5ad(
            args.results_dir / "mixed_6k_cellornucem.h5ad", compression="gzip"
        )

    print(json.dumps(summary, indent=2))
    print(f"\nResults: {args.results_dir.resolve()}")
    print(f"Figure:  {args.figure.resolve()}")


if __name__ == "__main__":
    main()
