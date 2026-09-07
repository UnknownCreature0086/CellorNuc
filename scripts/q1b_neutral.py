#!/usr/bin/env python3
"""Explain Neutral CellorNucEM calls using signature-count distributions.

This script reruns the Q1 global CellorNucEM analysis through the public API,
then assigns each Neutral droplet to the exact source-code condition that kept
it Neutral: too few signature counts, an intermediate posterior, or an
unrecognized assay label.
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


DEFAULT_INPUT = PROJECT_ROOT / "data/mixed_mBDRC_5pCells_multiomeNuclei_6k.h5ad"
DEFAULT_RESULTS = PROJECT_ROOT / "results/q1"
DEFAULT_FIGURE = PROJECT_ROOT / "figures/q1b_neutral_origins.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--min-counts", type=int, default=10)
    parser.add_argument("--gamma-hi", type=float, default=0.95)
    parser.add_argument("--gamma-lo", type=float, default=0.05)
    return parser.parse_args()


def describe_counts(obs: pd.DataFrame, min_counts: int) -> pd.DataFrame:
    """Return interpretable distribution statistics by known modality."""
    quantiles = {
        "q01": 0.01,
        "q05": 0.05,
        "q10": 0.10,
        "q25": 0.25,
        "median": 0.50,
        "q75": 0.75,
        "q90": 0.90,
        "q95": 0.95,
        "q99": 0.99,
    }
    rows: list[dict[str, float | int | str]] = []
    for source in ["cell", "nucleus", "all"]:
        values = (
            obs["modality_counts"]
            if source == "all"
            else obs.loc[obs["source"] == source, "modality_counts"]
        )
        below = values < min_counts
        row: dict[str, float | int | str] = {
            "source": source,
            "n": int(values.size),
            "mean": float(values.mean()),
            "sd": float(values.std()),
            "min": int(values.min()),
            "max": int(values.max()),
        }
        row.update({name: float(values.quantile(q)) for name, q in quantiles.items()})
        row["below_min_counts_n"] = int(below.sum())
        row["below_min_counts_percent"] = float(100 * below.mean())
        rows.append(row)
    return pd.DataFrame(rows).set_index("source")


def make_count_bins(obs: pd.DataFrame) -> pd.DataFrame:
    """Tabulate the low end of modality_counts without hiding zeros."""
    labels = ["0", "1-4", "5-9", "10-19", "20-49", "50-99", "100+"]
    count_bin = pd.cut(
        obs["modality_counts"],
        bins=[-1, 0, 4, 9, 19, 49, 99, np.inf],
        labels=labels,
    )
    counts = pd.crosstab(obs["source"], count_bin).reindex(
        index=["cell", "nucleus"], columns=labels, fill_value=0
    )
    long = counts.rename_axis(index="source", columns="count_bin").stack().rename(
        "n_droplets"
    ).reset_index()
    source_totals = obs["source"].value_counts()
    long["percent_within_source"] = [
        100 * row.n_droplets / source_totals[row.source]
        for row in long.itertuples(index=False)
    ]
    return long


def classify_neutral_causes(
    obs: pd.DataFrame,
    min_counts: int,
    gamma_lo: float,
    gamma_hi: float,
) -> pd.Series:
    """Mirror the mixed-dataset classification gates in CellorNucEM."""
    neutral = obs["classification"] == "Neutral"
    recognized = obs["source"].isin(["cell", "nucleus"])
    enough = obs["modality_counts"] >= min_counts
    intermediate = obs["p_cell"].between(gamma_lo, gamma_hi, inclusive="both")

    cause = pd.Series("Not Neutral", index=obs.index, dtype=object)
    cause.loc[neutral & ~recognized] = "Unrecognized assay label"
    cause.loc[neutral & recognized & ~enough] = "Below min_counts"
    cause.loc[neutral & recognized & enough & intermediate] = "Intermediate posterior"
    cause.loc[neutral & (cause == "Not Neutral")] = "Other/invalid state"
    return cause


def neutral_cause_table(obs: pd.DataFrame) -> pd.DataFrame:
    """Return cause counts and percentages for Neutral droplets."""
    neutral = obs.loc[obs["neutral_cause"] != "Not Neutral"].copy()
    sources = ["cell", "nucleus", "all"]
    causes = [
        "Below min_counts",
        "Intermediate posterior",
        "Unrecognized assay label",
        "Other/invalid state",
    ]
    rows: list[dict[str, float | int | str]] = []
    for source in sources:
        source_all = obs if source == "all" else obs.loc[obs["source"] == source]
        source_neutral = (
            neutral if source == "all" else neutral.loc[neutral["source"] == source]
        )
        for cause in causes:
            n = int((source_neutral["neutral_cause"] == cause).sum())
            rows.append(
                {
                    "source": source,
                    "neutral_cause": cause,
                    "n_droplets": n,
                    "percent_of_source": 100 * n / len(source_all),
                    "percent_of_neutral_in_source": (
                        100 * n / len(source_neutral) if len(source_neutral) else 0.0
                    ),
                }
            )
    return pd.DataFrame(rows)


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def plot_neutral_origins(
    obs: pd.DataFrame,
    output_path: Path,
    min_counts: int,
    gamma_lo: float,
    gamma_hi: float,
) -> None:
    """Visualize count distributions, the count gate, and Neutral causes."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    colors = {"cell": "#e76f51", "nucleus": "#457b9d"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)

    # Log1p keeps zero-count droplets visible while accommodating the long tail.
    transformed_bins = np.linspace(
        0, np.log10(obs["modality_counts"].max() + 1), 55
    )
    for source, color in colors.items():
        values = np.log10(obs.loc[obs["source"] == source, "modality_counts"] + 1)
        axes[0, 0].hist(
            values,
            bins=transformed_bins,
            density=True,
            histtype="step",
            linewidth=2,
            color=color,
            label=source,
        )
    axes[0, 0].axvline(
        np.log10(min_counts + 1), color="black", linestyle="--", linewidth=1.5
    )
    axes[0, 0].set(
        title="Distribution of signature evidence",
        xlabel="log10(modality_counts + 1)",
        ylabel="Density",
    )
    axes[0, 0].legend(frameon=False)

    for source, color in colors.items():
        values = obs.loc[obs["source"] == source, "modality_counts"].to_numpy()
        x, y = ecdf(values)
        axes[0, 1].step(x + 1, y, where="post", color=color, linewidth=2, label=source)
    axes[0, 1].axvline(
        min_counts + 1, color="black", linestyle="--", linewidth=1.5,
        label=f"min_counts = {min_counts}",
    )
    axes[0, 1].set_xscale("log")
    axes[0, 1].set(
        title="Empirical cumulative distribution",
        xlabel="modality_counts + 1 (log scale)",
        ylabel="Fraction at or below x",
        ylim=(0, 1.01),
    )
    axes[0, 1].legend(frameon=False)

    causes = ["Below min_counts", "Intermediate posterior"]
    bottom = np.zeros(2)
    x_positions = np.arange(2)
    cause_colors = ["#8d99ae", "#f4a261"]
    for cause, color in zip(causes, cause_colors):
        values = np.array(
            [
                ((obs["source"] == source) & (obs["neutral_cause"] == cause)).sum()
                for source in ["cell", "nucleus"]
            ]
        )
        axes[1, 0].bar(
            x_positions, values, bottom=bottom, color=color, label=cause
        )
        for x_position, base, value in zip(x_positions, bottom, values):
            if value:
                axes[1, 0].text(
                    x_position, base + value / 2, str(value),
                    ha="center", va="center", fontsize=10,
                )
        bottom += values
    axes[1, 0].set_xticks(x_positions, ["cell", "nucleus"])
    axes[1, 0].set(
        title="Why Neutral droplets remained Neutral",
        ylabel="Number of Neutral droplets",
    )
    axes[1, 0].legend(frameon=False)

    neutral = obs.loc[obs["classification"] == "Neutral"]
    for source, color in colors.items():
        group = neutral.loc[neutral["source"] == source]
        axes[1, 1].scatter(
            group["modality_counts"] + 1,
            group["p_cell"],
            s=10,
            alpha=0.5,
            linewidths=0,
            color=color,
            label=source,
        )
    axes[1, 1].axvspan(1, min_counts, color="#8d99ae", alpha=0.15)
    axes[1, 1].axvline(
        min_counts + 1, color="black", linestyle="--", linewidth=1.2
    )
    axes[1, 1].axhspan(gamma_lo, gamma_hi, color="#f4a261", alpha=0.12)
    axes[1, 1].axhline(gamma_lo, color="black", linestyle=":", linewidth=1)
    axes[1, 1].axhline(gamma_hi, color="black", linestyle=":", linewidth=1)
    axes[1, 1].set_xscale("log")
    axes[1, 1].set(
        title="Neutral calls under the two classification gates",
        xlabel="modality_counts + 1 (log scale)",
        ylabel="p_cell",
        ylim=(-0.02, 1.02),
    )
    axes[1, 1].legend(frameon=False)

    fig.suptitle("Q1b: origin of Neutral CellorNucEM calls", fontsize=15)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    sc_genes, sn_genes = Scanpyplus.LoadGeneSignatures()
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

    obs = pd.DataFrame(index=result.obs_names)
    obs["source"] = result.obs["suspension_type"].astype(str)
    obs["modality_counts"] = result.obs["modality_counts"].astype(int)
    obs["sc_frac"] = result.obs["sc_frac"].astype(float)
    obs["p_cell"] = result.obs["p_cell"].astype(float)
    obs["classification"] = result.obs["modality_classification"].astype(str)
    obs["neutral_cause"] = classify_neutral_causes(
        obs, args.min_counts, args.gamma_lo, args.gamma_hi
    )

    distribution = describe_counts(obs, args.min_counts)
    count_bins = make_count_bins(obs)
    causes = neutral_cause_table(obs)
    neutral = obs.loc[obs["classification"] == "Neutral"]

    cause_totals = neutral["neutral_cause"].value_counts()
    low_count_neutral = int(cause_totals.get("Below min_counts", 0))
    intermediate_neutral = int(cause_totals.get("Intermediate posterior", 0))
    summary = {
        "parameters": {
            "min_counts": args.min_counts,
            "gamma_lo": args.gamma_lo,
            "gamma_hi": args.gamma_hi,
            "fit": "global",
            "layer": "counts",
        },
        "definition": (
            "modality_counts is the per-droplet sum of raw counts over the "
            "sc and sn signature genes found in the dataset; it is not total UMIs."
        ),
        "n_droplets": int(len(obs)),
        "n_neutral": int(len(neutral)),
        "neutral_percent": float(100 * len(neutral) / len(obs)),
        "neutral_causes": {
            "below_min_counts": low_count_neutral,
            "below_min_counts_percent_of_neutral": float(
                100 * low_count_neutral / len(neutral)
            ),
            "intermediate_posterior": intermediate_neutral,
            "intermediate_posterior_percent_of_neutral": float(
                100 * intermediate_neutral / len(neutral)
            ),
            "unrecognized_assay_label": int(
                cause_totals.get("Unrecognized assay label", 0)
            ),
            "other_invalid_state": int(cause_totals.get("Other/invalid state", 0)),
        },
        "by_source": {
            source: {
                "median_modality_counts": float(distribution.loc[source, "median"]),
                "below_min_counts_n": int(
                    distribution.loc[source, "below_min_counts_n"]
                ),
                "below_min_counts_percent": float(
                    distribution.loc[source, "below_min_counts_percent"]
                ),
                "neutral_n": int((neutral["source"] == source).sum()),
                "neutral_below_min_counts_n": int(
                    (
                        (neutral["source"] == source)
                        & (neutral["neutral_cause"] == "Below min_counts")
                    ).sum()
                ),
                "neutral_intermediate_posterior_n": int(
                    (
                        (neutral["source"] == source)
                        & (neutral["neutral_cause"] == "Intermediate posterior")
                    ).sum()
                ),
            }
            for source in ["cell", "nucleus"]
        },
    }

    distribution.to_csv(args.results_dir / "q1b_modality_counts_summary.csv")
    count_bins.to_csv(args.results_dir / "q1b_modality_counts_bins.csv", index=False)
    causes.to_csv(args.results_dir / "q1b_neutral_causes.csv", index=False)
    neutral.to_csv(
        args.results_dir / "q1b_neutral_droplets.csv.gz", compression="gzip"
    )
    with (args.results_dir / "q1b_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    plot_neutral_origins(
        obs, args.figure, args.min_counts, args.gamma_lo, args.gamma_hi
    )

    print(json.dumps(summary, indent=2))
    print(f"\nTables: {args.results_dir.resolve()}")
    print(f"Figure: {args.figure.resolve()}")


if __name__ == "__main__":
    main()
