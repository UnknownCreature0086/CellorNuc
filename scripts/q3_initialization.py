#!/usr/bin/env python3
"""Compare anchored and split initialization for CellorNucEM on dataset B.

The two public API calls use the same observations, signatures, raw-count
layer, and all CellorNucEM model defaults. Only ``init`` changes. The script
also reconstructs the initialization values from the source helpers so that
the word "anchor" can be connected to the actual starting parameters.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path
from typing import Any

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
DEFAULT_RESULTS = PROJECT_ROOT / "results/q3"
DEFAULT_FIGURE = PROJECT_ROOT / "figures/q3_initialization_comparison.png"

CLASS_ORDER = ["Cell-like", "Nucleus-like", "Neutral"]
SCHEME_COLORS = {"anchored": "#2a9d8f", "split": "#7b2cbf"}
COMPONENT_COLORS = {"cell": "#e76f51", "nucleus": "#457b9d"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument(
        "--large-absolute",
        type=float,
        default=0.05,
        help="Practical large-difference cutoff on the 0-1 component-mean scale.",
    )
    parser.add_argument(
        "--large-relative",
        type=float,
        default=0.10,
        help="Large-difference cutoff relative to anchored component separation.",
    )
    parser.add_argument(
        "--convergence-audit-iterations",
        type=int,
        default=0,
        help=(
            "Optionally rerun the private mixture helper for this many iterations "
            "to test whether default-fit differences persist (0 disables)."
        ),
    )
    return parser.parse_args()


def component_statistics(ab: tuple[float, float]) -> dict[str, float]:
    """Return interpretable summaries of Beta(a,b)."""
    a, b = (float(ab[0]), float(ab[1]))
    kappa = a + b
    p = a / kappa
    return {
        "a": a,
        "b": b,
        "p": p,
        "kappa": kappa,
        "latent_sd": float(np.sqrt(p * (1 - p) / (kappa + 1))),
        "rho": 1 / (kappa + 1),
    }


def python_value(value: Any) -> Any:
    """Recursively convert NumPy values to JSON-safe Python values."""
    if isinstance(value, dict):
        return {str(key): python_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [python_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return [python_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def count_inputs(
    adata: sc.AnnData,
    sc_genes: list[str],
    sn_genes: list[str],
    layer: str,
    min_counts: int = 10,
) -> dict[str, np.ndarray]:
    """Create the exact compressed arrays used by CellorNucEM."""
    x = np.round(Scanpyplus._set_counts_EM(adata, sc_genes, layer)).astype(int)
    y = np.round(Scanpyplus._set_counts_EM(adata, sn_genes, layer)).astype(int)
    m = x + y
    enough = m >= min_counts
    pairs, counts = np.unique(
        np.stack([x[enough], m[enough]], axis=1), axis=0, return_counts=True
    )
    return {
        "x": x,
        "y": y,
        "m": m,
        "enough": enough,
        "xu": pairs[:, 0],
        "mu": pairs[:, 1],
        "counts": counts,
    }


def initialization_table(inputs: dict[str, np.ndarray]) -> tuple[pd.DataFrame, float]:
    """Reconstruct lines 1583-1603 of the supplied implementation."""
    xu, mu, counts = inputs["xu"], inputs["mu"], inputs["counts"]
    observed_p = xu / np.maximum(mu, 1)
    split_cutoff = float(np.average(observed_p, weights=counts * mu))
    n = int(counts.sum())
    rows: list[dict[str, Any]] = []

    schemes = {
        "anchored": {
            "cell_mask": observed_p > 0.7,
            "nucleus_mask": observed_p < 0.3,
            "cell_rule": "sc_frac > 0.7",
            "nucleus_rule": "sc_frac < 0.3",
        },
        "split": {
            "cell_mask": observed_p >= split_cutoff,
            "nucleus_mask": observed_p < split_cutoff,
            "cell_rule": f"sc_frac >= weighted mean ({split_cutoff:.6f})",
            "nucleus_rule": f"sc_frac < weighted mean ({split_cutoff:.6f})",
        },
    }

    for scheme, spec in schemes.items():
        high = spec["cell_mask"]
        low = spec["nucleus_mask"]
        if scheme == "anchored":
            initial_lambda = (
                float(np.clip(counts[high].sum() / n, 0.01, 0.99))
                if high.sum()
                else 0.9
            )
        else:
            initial_lambda = float(counts[high].sum() / n)
        for component, mask, rule, fallback in [
            ("cell", high, spec["cell_rule"], (19.0, 1.0)),
            ("nucleus", low, spec["nucleus_rule"], (1.0, 19.0)),
        ]:
            if scheme == "anchored" and int(mask.sum()) < 2:
                ab = fallback
                used_fallback = True
            else:
                ab = Scanpyplus._fit_bb_weighted_EM(
                    xu[mask], mu[mask], counts[mask]
                )
                used_fallback = False
            rows.append(
                {
                    "scheme": scheme,
                    "component": component,
                    "subset_rule": rule,
                    "subset_n_unique_pairs": int(mask.sum()),
                    "subset_n_droplets": int(counts[mask].sum()),
                    "initial_lambda_cell": initial_lambda,
                    "used_fallback": used_fallback,
                    **{
                        f"initial_{key}": value
                        for key, value in component_statistics(ab).items()
                    },
                }
            )
    return pd.DataFrame(rows), split_cutoff


def run_fit(
    adata: sc.AnnData,
    sc_genes: list[str],
    sn_genes: list[str],
    init: str,
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Run the public API in place and preserve only the comparison outputs."""
    # layer and print_mode are data-selection/reporting choices. Every
    # statistical model argument other than init is left at its API default.
    Scanpyplus.CellorNucEM(
        adata,
        sc_genes,
        sn_genes,
        layer="counts",
        init=init,
        print_mode=False,
        copy=False,
    )
    fit = dict(adata.uns["modality_em"]["global"])
    state = adata.obs[
        ["modality_counts", "sc_frac", "p_cell", "modality_classification"]
    ].copy()
    state["modality_classification"] = state["modality_classification"].astype(str)
    return fit, state


def fitted_parameter_table(fits: dict[str, dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scheme, fit in fits.items():
        for component, ab_key, weight in [
            ("cell", "ab_cell", float(fit["lam_cell"])),
            ("nucleus", "ab_nuc", 1 - float(fit["lam_cell"])),
        ]:
            rows.append(
                {
                    "scheme": scheme,
                    "component": component,
                    **component_statistics(fit[ab_key]),
                    "mixture_weight": weight,
                    "loglik": float(fit["loglik"]),
                    "bic_1": float(fit["bic_1"]),
                    "bic_2": float(fit["bic_2"]),
                    "delta_bic_1_minus_2": float(fit["bic_1"] - fit["bic_2"]),
                }
            )
    return pd.DataFrame(rows)


def convergence_audit_table(
    inputs: dict[str, np.ndarray], max_iter: int
) -> pd.DataFrame:
    """Run a longer private-helper fit and record its actual stopping iteration."""
    rows: list[dict[str, Any]] = []
    for scheme in ["anchored", "split"]:
        trace = io.StringIO()
        with contextlib.redirect_stdout(trace):
            fit = Scanpyplus._fit_two_comp_EM(
                inputs["xu"],
                inputs["mu"],
                inputs["counts"],
                init=scheme,
                max_iter=max_iter,
                verbose=True,
            )
        iteration_lines = [
            line for line in trace.getvalue().splitlines() if line.startswith("EM iter")
        ]
        last_iteration = (
            int(iteration_lines[-1].split()[2]) if iteration_lines else -1
        )
        cell = component_statistics(fit["ab_cell"])
        nucleus = component_statistics(fit["ab_nuc"])
        rows.append(
            {
                "scheme": scheme,
                "requested_max_iter": max_iter,
                "last_iteration_evaluated_zero_based": last_iteration,
                "n_iterations_evaluated": last_iteration + 1,
                "stopped_before_cap": (last_iteration + 1) < max_iter,
                "cell_p": cell["p"],
                "nucleus_p": nucleus["p"],
                "cell_kappa": cell["kappa"],
                "nucleus_kappa": nucleus["kappa"],
                "lambda_cell": float(fit["lam_cell"]),
                "loglik": float(fit["loglik"]),
            }
        )
    return pd.DataFrame(rows)


def mean_comparison_table(
    parameters: pd.DataFrame,
    absolute_cutoff: float,
    relative_cutoff: float,
) -> pd.DataFrame:
    p = parameters.pivot(index="component", columns="scheme", values="p")
    anchored = parameters.loc[parameters["scheme"] == "anchored"].set_index(
        "component"
    )
    anchored_separation = float(
        anchored.loc["cell", "p"] - anchored.loc["nucleus", "p"]
    )
    result = pd.DataFrame(
        {
            "component": ["cell", "nucleus"],
            "p_anchored": [p.loc["cell", "anchored"], p.loc["nucleus", "anchored"]],
            "p_split": [p.loc["cell", "split"], p.loc["nucleus", "split"]],
        }
    )
    result["signed_difference_split_minus_anchored"] = (
        result["p_split"] - result["p_anchored"]
    )
    result["absolute_difference"] = result[
        "signed_difference_split_minus_anchored"
    ].abs()
    result["percentage_point_difference"] = 100 * result["absolute_difference"]
    result["fraction_of_anchored_separation"] = (
        result["absolute_difference"] / anchored_separation
    )
    result["large_absolute_cutoff"] = absolute_cutoff
    result["large_relative_cutoff"] = relative_cutoff
    result["large_by_prespecified_rule"] = (
        (result["absolute_difference"] >= absolute_cutoff)
        & (result["fraction_of_anchored_separation"] >= relative_cutoff)
    )
    return result


def classification_tables(
    states: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for scheme, state in states.items():
        counts = state["modality_classification"].value_counts().reindex(
            CLASS_ORDER, fill_value=0
        )
        for label, count in counts.items():
            rows.append(
                {
                    "scheme": scheme,
                    "classification": label,
                    "n_droplets": int(count),
                    "fraction": float(count / len(state)),
                    "percent": float(100 * count / len(state)),
                }
            )
    count_table = pd.DataFrame(rows)
    crosstab = pd.crosstab(
        pd.Categorical(
            states["anchored"]["modality_classification"], categories=CLASS_ORDER
        ),
        pd.Categorical(
            states["split"]["modality_classification"], categories=CLASS_ORDER
        ),
        dropna=False,
    )
    crosstab.index.name = "anchored_classification"
    crosstab.columns.name = "split_classification"
    return count_table, crosstab


def posterior_summary(
    states: dict[str, pd.DataFrame], enough: np.ndarray
) -> pd.DataFrame:
    delta = (
        states["split"]["p_cell"].to_numpy()
        - states["anchored"]["p_cell"].to_numpy()
    )
    rows: list[dict[str, Any]] = []
    for population, mask in [
        ("all droplets", np.ones(len(delta), dtype=bool)),
        ("informative droplets", enough),
    ]:
        values = delta[mask]
        absolute = np.abs(values)
        rows.append(
            {
                "population": population,
                "n": int(mask.sum()),
                "mean_signed_difference": float(values.mean()),
                "mean_absolute_difference": float(absolute.mean()),
                "median_absolute_difference": float(np.median(absolute)),
                "q90_absolute_difference": float(np.quantile(absolute, 0.90)),
                "q95_absolute_difference": float(np.quantile(absolute, 0.95)),
                "max_absolute_difference": float(absolute.max()),
                "pearson_correlation": float(
                    np.corrcoef(
                        states["anchored"].loc[mask, "p_cell"],
                        states["split"].loc[mask, "p_cell"],
                    )[0, 1]
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_diagnostics(
    states: dict[str, pd.DataFrame],
    parameters: pd.DataFrame,
    crosstab: pd.DataFrame,
    enough: np.ndarray,
    split_cutoff: float,
    output_path: Path,
) -> None:
    """Create the Q3 data, parameter, posterior, and decision diagnostics."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), constrained_layout=True)

    ax = axes[0, 0]
    frac = states["anchored"].loc[enough, "sc_frac"].dropna()
    ax.hist(frac, bins=np.linspace(0, 1, 51), color="#a8adb3", edgecolor="white")
    ax.set_yscale("log")
    ax.axvspan(0, 0.3, color=COMPONENT_COLORS["nucleus"], alpha=0.12)
    ax.axvspan(0.7, 1, color=COMPONENT_COLORS["cell"], alpha=0.10)
    ax.axvline(
        split_cutoff,
        color=SCHEME_COLORS["split"],
        linestyle=":",
        linewidth=2,
        label=f"split cutoff = {split_cutoff:.3f}",
    )
    for _, row in parameters.iterrows():
        ax.axvline(
            row["p"],
            color=COMPONENT_COLORS[row["component"]],
            linestyle="-" if row["scheme"] == "anchored" else "--",
            linewidth=2,
            label=f"{row['scheme']} {row['component']} p = {row['p']:.3f}",
        )
    ax.set(
        title="A. Same sc-fraction data, different starting partitions",
        xlabel="sc_frac among informative droplets",
        ylabel="Number of droplets (log scale)",
        xlim=(0, 1),
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, frameon=False, fontsize=8, loc="upper left")

    ax = axes[0, 1]
    for component in ["cell", "nucleus"]:
        part = parameters.loc[parameters["component"] == component].set_index("scheme")
        values = [part.loc["anchored", "p"], part.loc["split", "p"]]
        ax.plot(
            [0, 1],
            values,
            marker="o",
            markersize=8,
            linewidth=2.5,
            color=COMPONENT_COLORS[component],
            label=f"{component} component",
        )
        for x_pos, value in enumerate(values):
            ax.text(x_pos, value + 0.025, f"{value:.3f}", ha="center", fontsize=9)
    ax.set_xticks([0, 1], ["anchored", "split"])
    ax.set_ylim(0, 1)
    ax.set(title="B. Final component means", ylabel="p = a / (a + b)")
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1, 0]
    hb = ax.hexbin(
        states["anchored"].loc[enough, "p_cell"],
        states["split"].loc[enough, "p_cell"],
        gridsize=45,
        bins="log",
        mincnt=1,
        cmap="viridis",
    )
    ax.plot([0, 1], [0, 1], color="white", linestyle="--", linewidth=1.5)
    ax.set(
        title="C. Posterior membership for the same droplets",
        xlabel="p_cell from anchored fit",
        ylabel="p_cell from split fit",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    fig.colorbar(hb, ax=ax, label="Bin count (log color scale)")

    ax = axes[1, 1]
    matrix = crosstab.reindex(index=CLASS_ORDER, columns=CLASS_ORDER, fill_value=0)
    image = ax.imshow(matrix.to_numpy(), cmap="Blues")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix.iloc[i, j])
            color = "white" if value > matrix.to_numpy().max() * 0.45 else "black"
            ax.text(j, i, f"{value:,}", ha="center", va="center", color=color)
    ax.set_xticks(range(len(CLASS_ORDER)), CLASS_ORDER, rotation=25, ha="right")
    ax.set_yticks(range(len(CLASS_ORDER)), CLASS_ORDER)
    ax.set(
        title="D. Hard-class agreement",
        xlabel="split initialization",
        ylabel="anchored initialization",
    )
    fig.colorbar(image, ax=ax, label="Number of droplets")

    fig.suptitle(
        "Q3: CellorNucEM initialization sensitivity on single-modality dataset B",
        fontsize=15,
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        raise ValueError("Dataset B must provide raw integer counts in layers['counts']")
    sc_genes, sn_genes = Scanpyplus.LoadGeneSignatures()
    inputs = count_inputs(adata, sc_genes, sn_genes, layer="counts")
    seeds, split_cutoff = initialization_table(inputs)

    fits: dict[str, dict[str, Any]] = {}
    states: dict[str, pd.DataFrame] = {}
    for scheme in ["anchored", "split"]:
        fits[scheme], states[scheme] = run_fit(
            adata, sc_genes, sn_genes, init=scheme
        )

    parameters = fitted_parameter_table(fits)
    mean_comparison = mean_comparison_table(
        parameters, args.large_absolute, args.large_relative
    )
    class_counts, class_crosstab = classification_tables(states)
    posterior = posterior_summary(states, inputs["enough"])
    convergence_audit = (
        convergence_audit_table(inputs, args.convergence_audit_iterations)
        if args.convergence_audit_iterations > 0
        else None
    )

    comparison = pd.DataFrame(index=adata.obs_names)
    comparison["modality_counts"] = states["anchored"]["modality_counts"]
    comparison["sc_frac"] = states["anchored"]["sc_frac"]
    comparison["informative"] = inputs["enough"]
    for scheme in ["anchored", "split"]:
        comparison[f"p_cell_{scheme}"] = states[scheme]["p_cell"]
        comparison[f"classification_{scheme}"] = states[scheme][
            "modality_classification"
        ]
    comparison["p_cell_difference_split_minus_anchored"] = (
        comparison["p_cell_split"] - comparison["p_cell_anchored"]
    )
    comparison["classification_changed"] = (
        comparison["classification_split"] != comparison["classification_anchored"]
    )

    seeds.to_csv(args.results_dir / "initialization_seeds.csv", index=False)
    parameters.to_csv(args.results_dir / "component_parameters.csv", index=False)
    mean_comparison.to_csv(
        args.results_dir / "component_mean_comparison.csv", index=False
    )
    class_counts.to_csv(args.results_dir / "classification_counts.csv", index=False)
    class_crosstab.to_csv(args.results_dir / "classification_crosstab.csv")
    posterior.to_csv(args.results_dir / "posterior_difference_summary.csv", index=False)
    if convergence_audit is not None:
        convergence_audit.to_csv(
            args.results_dir / "convergence_audit.csv", index=False
        )
    comparison.to_csv(
        args.results_dir / "droplet_comparison.csv.gz", compression="gzip"
    )

    plot_diagnostics(
        states,
        parameters,
        class_crosstab,
        inputs["enough"],
        split_cutoff,
        args.figure,
    )

    disagreement = comparison["classification_changed"].to_numpy()
    enough = inputs["enough"]
    summary = {
        "input": str(args.input.resolve()),
        "shape": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
        "comparison_design": {
            "only_model_argument_changed": "init",
            "schemes": ["anchored", "split"],
            "raw_count_layer": "counts",
            "other_CellorNucEM_model_parameters": "API defaults",
            "effective_defaults": {
                "min_counts": 10,
                "gamma_hi": 0.9,
                "gamma_lo": 0.1,
                "anchor_lo": 0.3,
                "anchor_hi": 0.7,
                "max_p_nuc": None,
                "min_p_cell": None,
            },
        },
        "signature_overlap": {
            "sc_found": [gene for gene in sc_genes if gene in adata.var_names],
            "sn_found": [gene for gene in sn_genes if gene in adata.var_names],
        },
        "data": {
            "n_informative": int(enough.sum()),
            "split_cutoff": split_cutoff,
        },
        "fits": python_value(fits),
        "mean_comparison": python_value(mean_comparison.to_dict(orient="records")),
        "impact": {
            "classification_disagreements_all": int(disagreement.sum()),
            "classification_disagreement_percent_all": float(100 * disagreement.mean()),
            "classification_disagreements_informative": int(
                (disagreement & enough).sum()
            ),
            "classification_disagreement_percent_informative": float(
                100 * disagreement[enough].mean()
            ),
            "posterior_difference": python_value(
                posterior.to_dict(orient="records")
            ),
        },
        "large_difference_rule": {
            "description": (
                "Flag a component when |delta p| is at least the absolute cutoff "
                "and at least the relative cutoff times the anchored separation."
            ),
            "absolute_cutoff": args.large_absolute,
            "relative_cutoff": args.large_relative,
            "is_a_formal_hypothesis_test": False,
        },
        "convergence_audit": (
            python_value(convergence_audit.to_dict(orient="records"))
            if convergence_audit is not None
            else None
        ),
        "outputs": {
            "figure_png": str(args.figure.resolve()),
            "figure_pdf": str(args.figure.with_suffix(".pdf").resolve()),
            "results_directory": str(args.results_dir.resolve()),
        },
    }
    with (args.results_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print("\nFinal component-mean comparison")
    print(mean_comparison.to_string(index=False))
    print("\nClassification counts")
    print(class_counts.to_string(index=False))
    print(
        f"\nChanged hard labels: {disagreement.sum():,}/{adata.n_obs:,} "
        f"({100 * disagreement.mean():.2f}%)"
    )
    if convergence_audit is not None:
        print("\nExtended-iteration convergence audit")
        print(convergence_audit.to_string(index=False))
    print(f"Saved tables to {args.results_dir.resolve()}")
    print(f"Saved figure to {args.figure.resolve()}")


if __name__ == "__main__":
    main()
