#!/usr/bin/env python3
"""Run CellorNucEM on the sc-only quiz dataset and summarize Q2.

The input has no assay column because every droplet came from an sc library.
Following the supplied demo, this script adds a constant ``modality='sc'``
column. CellorNucEM uses that column only to name the final classes; the EM
fit itself uses raw counts from the sc/sn signature genes.
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


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scanpyplus import Scanpyplus  # noqa: E402


DEFAULT_INPUT = PROJECT_ROOT / "data/single_KidneyRaji_sc_6k.h5ad"
DEFAULT_RESULTS = PROJECT_ROOT / "results/q2"
DEFAULT_HISTOGRAM = PROJECT_ROOT / "figures/q2_sc_frac_histogram.png"
DEFAULT_DECISION_FIGURE = PROJECT_ROOT / "figures/q2_score_to_classification.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--histogram", type=Path, default=DEFAULT_HISTOGRAM)
    parser.add_argument(
        "--decision-figure", type=Path, default=DEFAULT_DECISION_FIGURE
    )
    parser.add_argument("--min-counts", type=int, default=10)
    parser.add_argument("--gamma-hi", type=float, default=0.95)
    parser.add_argument("--gamma-lo", type=float, default=0.05)
    parser.add_argument("--max-p-nuc", type=float, default=0.2)
    parser.add_argument("--min-p-cell", type=float, default=0.8)
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


def make_h5ad_serializable(value: object) -> object:
    """Convert tuple-valued EM parameters into AnnData-writable arrays."""
    if isinstance(value, dict):
        return {key: make_h5ad_serializable(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return np.asarray(value)
    return value


def decision_category(
    obs: pd.DataFrame, min_counts: int, gamma_lo: float, gamma_hi: float
) -> pd.Categorical:
    """Expose the gates hidden by sc-only fallback-to-Typical labeling."""
    enough = obs["modality_counts"].to_numpy() >= min_counts
    p_cell = obs["p_cell"].to_numpy()
    values = np.full(len(obs), "Below min_counts", dtype=object)
    values[enough & (p_cell < gamma_lo)] = "Confident nucleus-like"
    values[enough & (p_cell > gamma_hi)] = "Confident cell-like"
    values[enough & (p_cell >= gamma_lo) & (p_cell <= gamma_hi)] = (
        "Intermediate posterior"
    )
    return pd.Categorical(
        values,
        categories=[
            "Confident nucleus-like",
            "Intermediate posterior",
            "Confident cell-like",
            "Below min_counts",
        ],
        ordered=True,
    )


def classification_summary(obs: pd.DataFrame) -> pd.DataFrame:
    order = ["Typical Cell (SC)", "Nucleus-like Cell (SC)"]
    counts = obs["modality_classification"].value_counts().reindex(order, fill_value=0)
    table = counts.rename("n_droplets").to_frame()
    table["fraction"] = table["n_droplets"] / len(obs)
    table["percent"] = 100 * table["fraction"]
    return table.rename_axis("modality_classification").reset_index()


def group_summary(obs: pd.DataFrame, group: str) -> pd.DataFrame:
    """Summarize whether nucleus-like calls concentrate in biological groups."""
    label = "Nucleus-like Cell (SC)"
    grouped = obs.groupby(group, observed=True)["modality_classification"]
    table = grouped.agg(
        n_droplets="size", nucleus_like_n=lambda x: int((x.astype(str) == label).sum())
    )
    table["nucleus_like_fraction"] = table["nucleus_like_n"] / table["n_droplets"]
    table["nucleus_like_percent"] = 100 * table["nucleus_like_fraction"]
    return table.sort_values(
        ["nucleus_like_percent", "n_droplets"], ascending=[False, False]
    ).reset_index()


def sc_frac_summary(obs: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    groups = [("All droplets", obs)] + list(
        obs.groupby("modality_classification", observed=True)
    )
    for label, group in groups:
        values = group["sc_frac"].dropna()
        rows.append(
            {
                "group": str(label),
                "n_droplets": int(len(group)),
                "n_defined_sc_frac": int(len(values)),
                "mean": float(values.mean()),
                "q05": float(values.quantile(0.05)),
                "q25": float(values.quantile(0.25)),
                "median": float(values.median()),
                "q75": float(values.quantile(0.75)),
                "q95": float(values.quantile(0.95)),
                "min": float(values.min()),
                "max": float(values.max()),
            }
        )
    return pd.DataFrame(rows)


def histogram_table(obs: pd.DataFrame, n_bins: int = 40) -> pd.DataFrame:
    """Save the exact counts underlying the requested histogram."""
    edges = np.linspace(0, 1, n_bins + 1)
    result = pd.DataFrame(
        {
            "bin_left": edges[:-1],
            "bin_right": edges[1:],
            "all_defined_sc_frac": np.histogram(obs["sc_frac"].dropna(), edges)[0],
        }
    )
    for label, column in [
        ("Typical Cell (SC)", "typical_cell"),
        ("Nucleus-like Cell (SC)", "nucleus_like_cell"),
    ]:
        values = obs.loc[obs["modality_classification"].astype(str) == label, "sc_frac"]
        result[column] = np.histogram(values.dropna(), edges)[0]
    return result


def plot_sc_frac_histogram(
    obs: pd.DataFrame,
    output_path: Path,
    cell_mean: float,
    nucleus_mean: float,
) -> None:
    """Plot the requested histogram using the same 40 bins as the demo."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    values = obs["sc_frac"].dropna()
    bins = np.linspace(0, 1, 41)
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.hist(values, bins=bins, color="#7b8fa1", edgecolor="white", linewidth=0.6)
    ax.axvline(
        nucleus_mean,
        color="#457b9d",
        linestyle="--",
        linewidth=1.8,
        label=f"Fitted nucleus mean = {nucleus_mean:.3f}",
    )
    ax.axvline(
        cell_mean,
        color="#e76f51",
        linestyle="--",
        linewidth=1.8,
        label=f"Fitted cell mean = {cell_mean:.3f}",
    )
    ax.set(
        title="Q2: sc-only dataset has a dominant cell-like mode\nand a smaller nucleus-like mode",
        xlabel="sc_frac = sc-signature counts / all signature counts",
        ylabel="Number of droplets",
        xlim=(0, 1),
    )
    ax.text(
        0.02,
        0.94,
        f"{len(values):,}/{len(obs):,} droplets have defined sc_frac",
        transform=ax.transAxes,
        va="top",
    )
    ax.legend(frameon=False)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_score_to_classification(
    obs: pd.DataFrame,
    output_path: Path,
    min_counts: int,
    gamma_lo: float,
    gamma_hi: float,
) -> None:
    """Connect raw sc_frac and posterior scores to the final hard classes."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    classification = obs["modality_classification"].astype(str)
    typical = classification == "Typical Cell (SC)"
    nucleus_like = classification == "Nucleus-like Cell (SC)"
    low_evidence = obs["modality_counts"].to_numpy() < min_counts

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)

    bins = np.linspace(0, 1, 41)
    axes[0].hist(
        [
            obs.loc[typical, "sc_frac"].dropna(),
            obs.loc[nucleus_like, "sc_frac"].dropna(),
        ],
        bins=bins,
        stacked=True,
        color=["#e76f51", "#457b9d"],
        label=["Typical Cell (SC)", "Nucleus-like Cell (SC)"],
        edgecolor="white",
        linewidth=0.35,
    )
    axes[0].set(
        title="Raw score distribution by hard classification",
        xlabel="sc_frac",
        ylabel="Number of droplets",
        xlim=(0, 1),
    )
    axes[0].legend(frameon=False)

    # Low-evidence droplets are drawn separately because their final Typical
    # label is an assay fallback, not a confident cell-like model decision.
    axes[1].scatter(
        obs.loc[low_evidence, "sc_frac"],
        obs.loc[low_evidence, "p_cell"],
        s=10,
        alpha=0.28,
        linewidths=0,
        color="#9e9e9e",
        label=f"Below min_counts ({low_evidence.sum():,})",
    )
    axes[1].scatter(
        obs.loc[~low_evidence & typical, "sc_frac"],
        obs.loc[~low_evidence & typical, "p_cell"],
        s=10,
        alpha=0.28,
        linewidths=0,
        color="#e76f51",
        label="Typical, enough counts",
    )
    axes[1].scatter(
        obs.loc[nucleus_like, "sc_frac"],
        obs.loc[nucleus_like, "p_cell"],
        s=13,
        alpha=0.65,
        linewidths=0,
        color="#457b9d",
        label="Nucleus-like hard call",
    )
    axes[1].axhline(gamma_lo, color="black", linestyle="--", linewidth=1)
    axes[1].axhline(gamma_hi, color="black", linestyle="--", linewidth=1)
    axes[1].set(
        title="Raw fraction, fitted posterior, and decision gates",
        xlabel="sc_frac",
        ylabel="p_cell",
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
    )
    axes[1].legend(frameon=False, loc="center left", fontsize=9)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        raise ValueError("Dataset has no raw-count layer named 'counts'.")
    if adata.n_obs == 0:
        raise ValueError("Dataset has no droplets.")

    # Dataset B is known to be sc-only. This metadata affects the displayed
    # class names but is explicitly not used by the EM likelihood.
    adata.obs["modality"] = "sc"
    sc_genes, sn_genes = Scanpyplus.LoadGeneSignatures()
    sc_found = [gene for gene in sc_genes if gene in adata.var_names]
    sn_found = [gene for gene in sn_genes if gene in adata.var_names]

    result = Scanpyplus.CellorNucEM(
        adata,
        sc_genes,
        sn_genes,
        assay_col="modality",
        sc_label="sc",
        sn_label="sn",
        layer="counts",
        min_counts=args.min_counts,
        gamma_hi=args.gamma_hi,
        gamma_lo=args.gamma_lo,
        max_p_nuc=args.max_p_nuc,
        min_p_cell=args.min_p_cell,
        print_mode=False,
        copy=True,
    )

    result.obs["decision_category"] = decision_category(
        result.obs, args.min_counts, args.gamma_lo, args.gamma_hi
    )
    classification = classification_summary(result.obs)
    decisions = (
        result.obs["decision_category"]
        .value_counts(sort=False)
        .rename("n_droplets")
        .to_frame()
    )
    decisions["fraction"] = decisions["n_droplets"] / result.n_obs
    decisions["percent"] = 100 * decisions["fraction"]
    decisions = decisions.rename_axis("decision_category").reset_index()

    label = "Nucleus-like Cell (SC)"
    nucleus_like_n = int(
        (result.obs["modality_classification"].astype(str) == label).sum()
    )
    enough = result.obs["modality_counts"].to_numpy() >= args.min_counts
    p_cell = result.obs["p_cell"].to_numpy()
    fit = result.uns["modality_em"]["global"]
    cell_mean = component_mean(fit["ab_cell"])
    nucleus_mean = component_mean(fit["ab_nuc"])
    expected_nucleus_informative = float(np.sum(1 - p_cell[enough]))

    hist_table = histogram_table(result.obs)
    frac_summary = sc_frac_summary(result.obs)
    by_cell_type = group_summary(result.obs, "reannotated_celltype")
    by_batch = group_summary(result.obs, "batch")

    classification.to_csv(args.results_dir / "classification_counts.csv", index=False)
    decisions.to_csv(args.results_dir / "decision_gate_counts.csv", index=False)
    hist_table.to_csv(args.results_dir / "sc_frac_histogram_bins.csv", index=False)
    frac_summary.to_csv(args.results_dir / "sc_frac_summary.csv", index=False)
    by_cell_type.to_csv(
        args.results_dir / "nucleus_like_by_cell_type.csv", index=False
    )
    by_batch.to_csv(args.results_dir / "nucleus_like_by_batch.csv", index=False)

    droplet_columns = [
        "batch",
        "leiden",
        "Predicted",
        "reannotated_celltype",
        "modality",
        "modality_counts",
        "sc_frac",
        "p_cell",
        "decision_category",
        "modality_fit",
        "modality_classification",
    ]
    result.obs[droplet_columns].to_csv(
        args.results_dir / "droplet_classifications.csv.gz", compression="gzip"
    )

    summary = {
        "input": str(args.input.resolve()),
        "shape": {"n_obs": int(result.n_obs), "n_vars": int(result.n_vars)},
        "parameters": {
            "fit": "global",
            "layer": "counts",
            "assay_col_added": "modality",
            "sc_label": "sc",
            "sn_label": "sn",
            "min_counts": args.min_counts,
            "gamma_hi": args.gamma_hi,
            "gamma_lo": args.gamma_lo,
            "init": "anchored",
            "max_p_nuc": args.max_p_nuc,
            "min_p_cell": args.min_p_cell,
        },
        "signature_overlap": {
            "sc_found": sc_found,
            "sc_missing": sorted(set(sc_genes).difference(sc_found)),
            "sn_found": sn_found,
            "sn_missing": sorted(set(sn_genes).difference(sn_found)),
        },
        "answer": {
            "nucleus_like_cell_n": nucleus_like_n,
            "nucleus_like_cell_fraction": nucleus_like_n / result.n_obs,
            "nucleus_like_cell_percent": 100 * nucleus_like_n / result.n_obs,
            "sc_frac_shape": "bimodal",
        },
        "evidence": {
            "informative_droplets": int(enough.sum()),
            "low_signature_count_droplets": int((~enough).sum()),
            "zero_signature_count_droplets": int(
                (result.obs["modality_counts"] == 0).sum()
            ),
            "defined_sc_frac_droplets": int(result.obs["sc_frac"].notna().sum()),
            "confident_nucleus_like": int((enough & (p_cell < args.gamma_lo)).sum()),
            "confident_cell_like": int((enough & (p_cell > args.gamma_hi)).sum()),
            "intermediate_posterior": int(
                (enough & (p_cell >= args.gamma_lo) & (p_cell <= args.gamma_hi)).sum()
            ),
        },
        "fit_diagnostics": {
            "cell_component_mean_sc_fraction": cell_mean,
            "nucleus_component_mean_sc_fraction": nucleus_mean,
            "component_separation": cell_mean - nucleus_mean,
            "lambda_cell": float(fit["lam_cell"]),
            "lambda_nucleus": float(1 - fit["lam_cell"]),
            "expected_nucleus_membership_informative": expected_nucleus_informative,
            "expected_nucleus_fraction_informative": (
                expected_nucleus_informative / int(enough.sum())
            ),
            "bic_1_component": float(fit["bic_1"]),
            "bic_2_component": float(fit["bic_2"]),
            "delta_bic_1_minus_2": float(fit["bic_1"] - fit["bic_2"]),
            "nucleus_mean_at_constraint": bool(
                np.isclose(nucleus_mean, args.max_p_nuc)
            ),
        },
        "outputs": {
            "histogram": str(args.histogram.resolve()),
            "decision_figure": str(args.decision_figure.resolve()),
        },
    }
    with (args.results_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    plot_sc_frac_histogram(
        result.obs, args.histogram, cell_mean=cell_mean, nucleus_mean=nucleus_mean
    )
    plot_score_to_classification(
        result.obs,
        args.decision_figure,
        min_counts=args.min_counts,
        gamma_lo=args.gamma_lo,
        gamma_hi=args.gamma_hi,
    )

    if not args.skip_h5ad:
        result.uns["modality_em"] = make_h5ad_serializable(
            result.uns["modality_em"]
        )
        result.write_h5ad(args.results_dir / "single_6k_cellornucem.h5ad")

    print(classification.to_string(index=False))
    print(
        f"\nQ2a: {nucleus_like_n:,}/{result.n_obs:,} = "
        f"{100 * nucleus_like_n / result.n_obs:.2f}% Nucleus-like Cell (SC)"
    )
    print("Q2b: sc_frac is bimodal (dominant high mode, smaller near-zero mode).")
    print(f"Saved results to {args.results_dir.resolve()}")
    print(f"Saved figures to {args.histogram.resolve()} and {args.decision_figure.resolve()}")


if __name__ == "__main__":
    main()
