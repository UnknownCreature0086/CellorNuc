#!/usr/bin/env python3
"""Run the Q3 initialization comparison on cell-only and nucleus-only dataset A."""

from __future__ import annotations

import argparse
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

from q3_initialization import (  # noqa: E402
    CLASS_ORDER,
    COMPONENT_COLORS,
    SCHEME_COLORS,
    classification_tables,
    component_statistics,
    convergence_audit_table,
    count_inputs,
    fitted_parameter_table,
    initialization_table,
    mean_comparison_table,
    posterior_summary,
    python_value,
    run_fit,
)
from scanpyplus import Scanpyplus  # noqa: E402


DEFAULT_INPUT = PROJECT_ROOT / "data/mixed_mBDRC_5pCells_multiomeNuclei_6k.h5ad"
DEFAULT_RESULTS = PROJECT_ROOT / "results/q3b"
DEFAULT_FIGURE = PROJECT_ROOT / "figures/q3b_dataset_a_single_modality.png"
DEFAULT_DATASET_B_SUMMARY = PROJECT_ROOT / "results/q3/summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument(
        "--dataset-b-summary", type=Path, default=DEFAULT_DATASET_B_SUMMARY
    )
    parser.add_argument("--large-absolute", type=float, default=0.05)
    parser.add_argument("--large-relative", type=float, default=0.10)
    parser.add_argument(
        "--convergence-audit-iterations",
        type=int,
        default=0,
        help="Optional long-run audit using the private EM helper (0 disables).",
    )
    return parser.parse_args()


def plot_subset_row(
    axes: np.ndarray,
    subset_name: str,
    states: dict[str, pd.DataFrame],
    parameters: pd.DataFrame,
    enough: np.ndarray,
    split_cutoff: float,
) -> None:
    """Plot raw evidence, final component means, and posterior sensitivity."""
    ax = axes[0]
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
        label=f"split cutoff={split_cutoff:.3f}",
    )
    for _, row in parameters.iterrows():
        ax.axvline(
            row["p"],
            color=COMPONENT_COLORS[row["component"]],
            linestyle="-" if row["scheme"] == "anchored" else "--",
            linewidth=1.8,
            label=f"{row['scheme']} {row['component']}={row['p']:.3f}",
        )
    ax.set(
        title=f"{subset_name}-only: sc_frac and fitted means",
        xlabel="sc_frac among informative droplets",
        ylabel="Droplets (log scale)",
        xlim=(0, 1),
    )
    ax.legend(frameon=False, fontsize=7, loc="upper left")

    ax = axes[1]
    for component in ["cell", "nucleus"]:
        part = parameters.loc[parameters["component"] == component].set_index("scheme")
        values = [part.loc["anchored", "p"], part.loc["split", "p"]]
        ax.plot(
            [0, 1],
            values,
            marker="o",
            markersize=7,
            linewidth=2.5,
            color=COMPONENT_COLORS[component],
            label=f"{component} component",
        )
        for x_pos, value in enumerate(values):
            ax.text(x_pos, value + 0.025, f"{value:.3f}", ha="center", fontsize=8)
    ax.set_xticks([0, 1], ["anchored", "split"])
    ax.set_ylim(0, 1)
    ax.set(title=f"{subset_name}-only: component means", ylabel="p = a/(a+b)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2]
    hb = ax.hexbin(
        states["anchored"].loc[enough, "p_cell"],
        states["split"].loc[enough, "p_cell"],
        gridsize=35,
        bins="log",
        mincnt=1,
        cmap="viridis",
    )
    ax.plot([0, 1], [0, 1], color="white", linestyle="--", linewidth=1.2)
    changed = (
        states["anchored"]["modality_classification"].to_numpy()
        != states["split"]["modality_classification"].to_numpy()
    )
    ax.text(
        0.03,
        0.97,
        f"hard labels changed: {changed.sum():,}/{len(changed):,}",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    ax.set(
        title=f"{subset_name}-only: posterior comparison",
        xlabel="p_cell, anchored",
        ylabel="p_cell, split",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    plt.colorbar(hb, ax=ax, label="Bin count (log color scale)")


def cross_dataset_table(
    subset_results: dict[str, dict[str, Any]],
    balanced_fits: dict[str, dict[str, Any]],
    dataset_b_summary: Path,
) -> pd.DataFrame:
    """Place dataset A single-source results beside the prior dataset B run."""
    records: list[dict[str, Any]] = []

    def append_fit(dataset: str, subset: str, fits: dict[str, dict[str, Any]]) -> None:
        for component, key in [("cell", "ab_cell"), ("nucleus", "ab_nuc")]:
            anchored = component_statistics(fits["anchored"][key])["p"]
            split = component_statistics(fits["split"][key])["p"]
            records.append(
                {
                    "dataset": dataset,
                    "single_modality_subset": subset,
                    "component": component,
                    "p_anchored": anchored,
                    "p_split": split,
                    "absolute_difference": abs(split - anchored),
                }
            )

    for subset, result in subset_results.items():
        append_fit("A (6k source subset)", subset, result["fits"])
    append_fit("A (balanced 6k control)", "cell+nucleus", balanced_fits)

    if dataset_b_summary.exists():
        with dataset_b_summary.open() as handle:
            dataset_b = json.load(handle)
        append_fit("B (6k)", "cell", dataset_b["fits"])
    return pd.DataFrame(records)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    if "suspension_type" not in adata.obs:
        raise ValueError("Dataset A must contain obs['suspension_type']")
    if "counts" not in adata.layers:
        raise ValueError("Dataset A must contain raw layers['counts']")
    observed_sources = set(adata.obs["suspension_type"].astype(str))
    if not {"cell", "nucleus"}.issubset(observed_sources):
        raise ValueError("suspension_type must contain both 'cell' and 'nucleus'")

    sc_genes, sn_genes = Scanpyplus.LoadGeneSignatures()
    subset_results: dict[str, dict[str, Any]] = {}
    combined_parameters: list[pd.DataFrame] = []
    combined_means: list[pd.DataFrame] = []
    combined_seeds: list[pd.DataFrame] = []
    combined_class_counts: list[pd.DataFrame] = []
    combined_transitions: list[pd.DataFrame] = []
    combined_posteriors: list[pd.DataFrame] = []
    combined_audits: list[pd.DataFrame] = []
    combined_droplets: list[pd.DataFrame] = []

    for subset_name in ["cell", "nucleus"]:
        mask = adata.obs["suspension_type"].astype(str) == subset_name
        subset = adata[mask].copy()
        inputs = count_inputs(subset, sc_genes, sn_genes, layer="counts")
        seeds, split_cutoff = initialization_table(inputs)
        fits: dict[str, dict[str, Any]] = {}
        states: dict[str, pd.DataFrame] = {}
        for scheme in ["anchored", "split"]:
            fits[scheme], states[scheme] = run_fit(
                subset, sc_genes, sn_genes, init=scheme
            )

        parameters = fitted_parameter_table(fits)
        means = mean_comparison_table(
            parameters, args.large_absolute, args.large_relative
        )
        class_counts, crosstab = classification_tables(states)
        posteriors = posterior_summary(states, inputs["enough"])
        audit = (
            convergence_audit_table(inputs, args.convergence_audit_iterations)
            if args.convergence_audit_iterations > 0
            else None
        )

        for table in [parameters, means, seeds, class_counts, posteriors]:
            table.insert(0, "single_modality_subset", subset_name)
        transition = crosstab.stack(future_stack=True).rename("n_droplets").reset_index()
        transition.insert(0, "single_modality_subset", subset_name)
        if audit is not None:
            audit.insert(0, "single_modality_subset", subset_name)

        droplets = pd.DataFrame(index=subset.obs_names)
        droplets["single_modality_subset"] = subset_name
        droplets["modality_counts"] = states["anchored"]["modality_counts"]
        droplets["sc_frac"] = states["anchored"]["sc_frac"]
        droplets["informative"] = inputs["enough"]
        for scheme in ["anchored", "split"]:
            droplets[f"p_cell_{scheme}"] = states[scheme]["p_cell"]
            droplets[f"classification_{scheme}"] = states[scheme][
                "modality_classification"
            ]
        droplets["classification_changed"] = (
            droplets["classification_anchored"] != droplets["classification_split"]
        )

        combined_parameters.append(parameters)
        combined_means.append(means)
        combined_seeds.append(seeds)
        combined_class_counts.append(class_counts)
        combined_transitions.append(transition)
        combined_posteriors.append(posteriors)
        if audit is not None:
            combined_audits.append(audit)
        combined_droplets.append(droplets)
        subset_results[subset_name] = {
            "n_obs": int(subset.n_obs),
            "n_informative": int(inputs["enough"].sum()),
            "split_cutoff": split_cutoff,
            "fits": fits,
            "states": states,
            "inputs": inputs,
            "parameters": parameters,
            "means": means,
            "classification_crosstab": crosstab,
            "posterior_summary": posteriors,
            "convergence_audit": audit,
        }

    parameters_all = pd.concat(combined_parameters, ignore_index=True)
    means_all = pd.concat(combined_means, ignore_index=True)
    seeds_all = pd.concat(combined_seeds, ignore_index=True)
    class_counts_all = pd.concat(combined_class_counts, ignore_index=True)
    transitions_all = pd.concat(combined_transitions, ignore_index=True)
    posteriors_all = pd.concat(combined_posteriors, ignore_index=True)
    droplets_all = pd.concat(combined_droplets)
    balanced_fits: dict[str, dict[str, Any]] = {}
    for scheme in ["anchored", "split"]:
        balanced_fits[scheme], _ = run_fit(adata, sc_genes, sn_genes, init=scheme)
    cross_dataset = cross_dataset_table(
        subset_results, balanced_fits, args.dataset_b_summary
    )

    parameters_all.to_csv(args.results_dir / "component_parameters.csv", index=False)
    means_all.to_csv(args.results_dir / "component_mean_comparison.csv", index=False)
    seeds_all.to_csv(args.results_dir / "initialization_seeds.csv", index=False)
    class_counts_all.to_csv(args.results_dir / "classification_counts.csv", index=False)
    transitions_all.to_csv(
        args.results_dir / "classification_transitions.csv", index=False
    )
    posteriors_all.to_csv(
        args.results_dir / "posterior_difference_summary.csv", index=False
    )
    droplets_all.to_csv(
        args.results_dir / "droplet_comparison.csv.gz", compression="gzip"
    )
    cross_dataset.to_csv(
        args.results_dir / "comparison_with_dataset_b.csv", index=False
    )
    if combined_audits:
        pd.concat(combined_audits, ignore_index=True).to_csv(
            args.results_dir / "convergence_audit.csv", index=False
        )

    fig, axes = plt.subplots(2, 3, figsize=(17, 10.5), constrained_layout=True)
    for row, subset_name in enumerate(["cell", "nucleus"]):
        result = subset_results[subset_name]
        plot_subset_row(
            axes[row],
            subset_name.capitalize(),
            result["states"],
            result["parameters"],
            result["inputs"]["enough"],
            result["split_cutoff"],
        )
    fig.suptitle(
        "Q3b: Initialization sensitivity after making dataset A single-modality",
        fontsize=15,
    )
    fig.savefig(args.figure, dpi=220, bbox_inches="tight")
    fig.savefig(args.figure.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    serializable_results = {}
    for subset_name, result in subset_results.items():
        serializable_results[subset_name] = {
            "n_obs": result["n_obs"],
            "n_informative": result["n_informative"],
            "split_cutoff": result["split_cutoff"],
            "fits": python_value(result["fits"]),
            "mean_comparison": python_value(
                result["means"].to_dict(orient="records")
            ),
            "classification_crosstab": python_value(
                result["classification_crosstab"].to_dict()
            ),
            "posterior_summary": python_value(
                result["posterior_summary"].to_dict(orient="records")
            ),
            "convergence_audit": (
                python_value(result["convergence_audit"].to_dict(orient="records"))
                if result["convergence_audit"] is not None
                else None
            ),
        }
    summary = {
        "input": str(args.input.resolve()),
        "shape": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
        "design": {
            "subsets": ["cell", "nucleus"],
            "only_model_argument_changed": "init",
            "schemes": ["anchored", "split"],
            "layer": "counts",
            "other_CellorNucEM_model_parameters": "API defaults",
            "large_difference_rule": {
                "absolute_cutoff": args.large_absolute,
                "relative_cutoff": args.large_relative,
            },
        },
        "signature_overlap": {
            "sc_found": [gene for gene in sc_genes if gene in adata.var_names],
            "sn_found": [gene for gene in sn_genes if gene in adata.var_names],
        },
        "subsets": serializable_results,
        "balanced_dataset_a_control": {
            "fits": python_value(balanced_fits),
            "component_means": {
                scheme: {
                    "cell": component_statistics(fit["ab_cell"])["p"],
                    "nucleus": component_statistics(fit["ab_nuc"])["p"],
                }
                for scheme, fit in balanced_fits.items()
            },
        },
        "outputs": {
            "results_directory": str(args.results_dir.resolve()),
            "figure_png": str(args.figure.resolve()),
            "figure_pdf": str(args.figure.with_suffix(".pdf").resolve()),
        },
    }
    with (args.results_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print("\nDefault-fit component means")
    print(
        means_all[
            [
                "single_modality_subset",
                "component",
                "p_anchored",
                "p_split",
                "absolute_difference",
                "large_by_prespecified_rule",
            ]
        ].to_string(index=False)
    )
    print("\nHard-class counts")
    print(class_counts_all.to_string(index=False))
    if combined_audits:
        print("\nExtended convergence audit")
        print(pd.concat(combined_audits, ignore_index=True).to_string(index=False))
    print(f"\nSaved tables to {args.results_dir.resolve()}")
    print(f"Saved figure to {args.figure.resolve()}")


if __name__ == "__main__":
    main()
