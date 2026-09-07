# CellorNuc starter quiz

Reproducible analyses for the CellorNuc starter quiz. Q1a evaluates the
CellorNucEM signature on an independent mixed cell/nucleus dataset; Q1b traces
the resulting Neutral calls to the count and posterior gates in the supplied
implementation. The datasets live in `data/`, and the professor-provided
ScanpyPlus toolkit lives in `scanpyplus/`.

## Environment

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate cellornuc
python -m ipykernel install --user --name cellornuc --display-name "Python (cellornuc)"
```

The local `.venv` is an equivalent ready-to-use environment when present:

```bash
source .venv/bin/activate
```

## Project layout

Generated files under `figures/` and `results/` are ignored by Git and can be
recreated from the scripts. Dataset files are also ignored because of their
size.

```text
CellorNuc/
├── README.md
│   Project overview, environment setup, file map, and run commands.
├── environment.yml
│   Conda environment for Scanpy, AnnData, scientific Python, and Jupyter.
├── URAP_quiz.pdf
│   Original quiz prompt/reference document.
├── .gitignore
│   Excludes environments, caches, datasets, and generated results/figures.
│
├── data/                                      # supplied locally; Git-ignored
│   ├── mixed_mBDRC_5pCells_multiomeNuclei_6k.h5ad
│   │   Q1 development dataset A: 3,000 cell + 3,000 nucleus droplets.
│   ├── mixed_mBDRC_5pCells_multiomeNuclei_full.h5ad
│   │   Full mixed dataset A, reserved for later confirmation.
│   ├── single_KidneyRaji_sc_6k.h5ad
│   │   Downsampled single-modality dataset B for Q2 development.
│   └── single_KidneyRaji_sc_full.h5ad
│       Full single-modality dataset B, reserved for later confirmation.
│
├── scanpyplus/                                # supplied analysis toolkit
│   ├── Scanpyplus.py
│   │   Professor-provided source, including LoadGeneSignatures and CellorNucEM.
│   ├── gene_signatures.json
│   │   Versioned sc/sn signature genes and their provenance.
│   └── CellorNuc_demo.ipynb
│       API examples used to align the quiz analyses with the intended workflow.
│
├── scripts/
│   ├── q1_mixed.py
│   │   Q1a: runs global CellorNucEM on dataset A; evaluates predictions,
│   │   coverage, AUROC, fit separation, and writes the annotated result.
│   ├── q1b_neutral.py
│   │   Q1b: reruns the same API call; summarizes modality_counts and assigns
│   │   every Neutral call to its exact count/posterior/metadata gate.
│   ├── inspect_data.py
│   │   Scaffold for reusable AnnData schema and count-layer inspection.
│   ├── q2_single.py
│   │   Scaffold for the single-modality Q2 analysis.
│   ├── q3_initialization.py
│   │   Scaffold for the Q3 initialization analysis.
│   └── bonus_celltypes.py
│       Scaffold for the optional cell-type-stratified analysis.
│
├── notebook/
│   └── CellorNuc_quiz.ipynb
│       Scaffold for assembling the final quiz narrative and executable results.
│
├── notes/
│   ├── source_reading.md
│   │   Detailed reading of CellorNucEM, its beta-binomial model, and API rules.
│   └── experiment_log.md
│       Submission-ready Q1a/Q1b answers, figure interpretations, parameters,
│       caveats, commands, and recorded findings.
│
├── figures/                                   # generated; Git-ignored
│   ├── .gitkeep
│   │   Keeps the otherwise-empty output directory in Git.
│   ├── q1_cellornucem_diagnostics.png
│   │   Q1a: known source vs prediction UMAPs plus sc_frac/p_cell distributions.
│   └── q1b_neutral_origins.png
│       Q1b: modality_counts histogram/ECDF, Neutral causes, and decision gates.
│
└── results/                                   # generated; Git-ignored
    ├── q1/
    │   ├── .gitkeep
    │   │   Keeps the Q1 result directory in Git.
    │   ├── summary.json
    │   │   Q1a parameters, signature overlap, fit diagnostics, and metrics.
    │   ├── classification_counts.csv
    │   │   Counts for the five metadata-aware CellorNucEM labels.
    │   ├── assay_by_classification.csv
    │   │   Known source crossed with metadata-aware CellorNucEM labels.
    │   ├── assay_by_prediction.csv
    │   │   Known source crossed with Cell-like/Nucleus-like/Neutral calls.
    │   ├── droplet_classifications.csv.gz
    │   │   Per-droplet Q1a evidence, posterior, prediction, and API label.
    │   ├── mixed_6k_cellornucem.h5ad
    │   │   Dataset A with CellorNucEM observation columns and fitted parameters.
    │   ├── q1b_summary.json
    │   │   Q1b headline distribution statistics and Neutral-cause totals.
    │   ├── q1b_modality_counts_summary.csv
    │   │   Count-distribution quantiles and below-threshold rates by source.
    │   ├── q1b_modality_counts_bins.csv
    │   │   Fine-grained modality_counts bins and percentages by source.
    │   ├── q1b_neutral_causes.csv
    │   │   Exact Neutral-cause counts and percentages by source.
    │   └── q1b_neutral_droplets.csv.gz
    │       Per-droplet records for all Neutral calls and their assigned cause.
    ├── q2/.gitkeep
    │   Placeholder for generated Q2 results.
    └── q3/.gitkeep
        Placeholder for generated Q3 results.
```

## Run Q1

Run commands from the repository root after activating the environment:

```bash
python scripts/q1_mixed.py
python scripts/q1b_neutral.py
```

`q1_mixed.py` runs the shared CellorNucEM fit and produces the Q1a evaluation.
`q1b_neutral.py` deliberately reruns the same public API call so its analysis is
independently reproducible, then explains every Neutral result from the API's
`min_counts` and posterior rules.

The current 6k run produced 4,783 correct calls among 4,849 confident calls
(98.64% accuracy at 80.82% overall coverage). Of 1,151 Neutral calls, 987
(85.75%) had fewer than 10 signature counts and 164 had an intermediate
posterior. See `notes/experiment_log.md` for the full interpretation and
submission-ready text.

Start with the 6k datasets for development; reserve the full datasets for the
later confirmation phase.
