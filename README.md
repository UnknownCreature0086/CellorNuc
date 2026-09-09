# CellorNuc Starter Quiz

This repository contains my analysis for the CellorNuc starter quiz.

## Start here

**Primary submission:** [`notebook/CellorNuc_quiz.ipynb`](notebook/CellorNuc_quiz.ipynb)

This is the concise, full-dataset notebook intended for review. It contains:

- a short mathematical and source-code overview of `CellorNucEM`;
- Q1 validation and Neutral-call analysis on full mixed dataset A;
- Q2 classification and the `sc_frac` histogram on full sc-only dataset B;
- Q3 anchored-versus-split initialization comparisons;
- a short bonus analysis of cell identity and classifier behavior; and
- the final conclusion and reflection.

The PDF attached to the submission email is a static export of this same
notebook. The `.ipynb` is included so that the code, saved outputs, and cell
ordering can also be inspected.

## Reading guide

| If you want to... | Open... |
|---|---|
| Read the submitted analysis | [`notebook/CellorNuc_quiz.ipynb`](notebook/CellorNuc_quiz.ipynb) |
| See the smaller development workflow (quite messy though) | [`notebook/CellorNuc_quiz_6k.ipynb`](notebook/CellorNuc_quiz_6k.ipynb) |
| Inspect the supplied implementation | [`scanpyplus/Scanpyplus.py`](scanpyplus/Scanpyplus.py) |
| See the original API example | [`scanpyplus/CellorNuc_demo.ipynb`](scanpyplus/CellorNuc_demo.ipynb) |
| Read the original quiz prompt | [`URAP_quiz.pdf`](URAP_quiz.pdf) |
| Inspect extended working notes | [`notes/source_reading.md`](notes/source_reading.md) and [`notes/experiment_log.md`](notes/experiment_log.md) |

The two files under `notes/` are AI-assisted working notes used to understand
the source and organize experiments. They are deliberately more extensive than
the submitted notebook and are not intended as separate quiz answers.

## Analysis workflow

The workflow was developed first on the two 6k datasets, where iteration and
debugging were faster. The same analysis was then simplified and run on the
full datasets for the primary notebook.

```text
6k datasets
    -> verify metadata, raw-count layer, signatures, and decision gates
    -> develop Q1-Q3 tables and figures
    -> inspect CellorNucEM mathematics and source behavior
    -> run and simplify the full-dataset notebook
    -> add a small cell-identity bonus analysis
```

## Repository structure

```text
CellorNuc/
├── README.md
├── URAP_quiz.pdf                     # original quiz prompt
├── environment.yml                  # reproducible Python environment
│
├── notebook/
│   ├── CellorNuc_quiz.ipynb         # primary full-dataset submission
│   └── CellorNuc_quiz_6k.ipynb      # development-scale notebook
│
├── scanpyplus/
│   ├── Scanpyplus.py                # supplied CellorNucEM implementation
│   ├── gene_signatures.json         # supplied sc/sn signatures
│   └── CellorNuc_demo.ipynb         # supplied/adapted API demonstration
│
├── scripts/                         # reproducible diagnostic analyses
│   ├── q1_mixed.py
│   ├── q1b_neutral.py
│   ├── q2_single.py
│   ├── q3_initialization.py
│   └── q3b_dataset_a_subsets.py
│
├── notes/                           # optional AI-assisted working notes
│   ├── source_reading.md
│   └── experiment_log.md
│
├── data/                            # supplied datasets; Git-ignored
├── figures/                         # generated diagnostics; Git-ignored
└── results/                         # generated tables/objects; Git-ignored
```

The data and generated outputs are excluded from Git because of their size.
The primary notebook retains the result tables and figures needed to read the
submission without opening those directories.

## Environment and reproduction

Create the environment:

```bash
conda env create -f environment.yml
conda activate cellornuc
```

Place the four supplied `.h5ad` files under `data/`, then open the primary
notebook from the repository root:

```bash
jupyter lab notebook/CellorNuc_quiz.ipynb
```

The notebook uses raw counts from `layers["counts"]` and contains saved outputs
for review. The scripts provide additional diagnostics but are not required to
follow the submitted answers.
