# CellorNuc starter quiz

Repository scaffold for the CellorNuc starter quiz. The supplied datasets live
in `data/`, and the professor-provided ScanpyPlus toolkit lives in
`scanpyplus/`.

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

- `scripts/`: one entry point per quiz component, plus data inspection
- `figures/`: generated plots
- `results/`: generated outputs grouped by question
- `notes/`: source-reading notes and experiment log
- `notebook/`: final quiz notebook

Start with the 6k datasets for development; reserve the full datasets for the
later confirmation phase.

