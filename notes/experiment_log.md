# Experiment log

Record each run's date, dataset, parameters, outputs, and observations here.

## 2026-09-07 — Q1a mixed dataset A (6k)

### Submission-ready answer

I ran a global `CellorNucEM` fit on the raw `counts` layer using the supplied
sc/sn signatures, `min_counts=10`, and posterior cutoffs of 0.95 and 0.05.
The dataset's `suspension_type` was passed as `assay_col`; the EM mixture does
not use this column when fitting, only afterward when naming agreement or
disagreement with the known assay source.

| Known source (`suspension_type`) | Cell-like call | Nucleus-like call | Neutral | Total |
|---|---:|---:|---:|---:|
| Cell | 2,735 | 59 | 206 | 3,000 |
| Nucleus | 7 | 2,048 | 945 | 3,000 |

No cell droplet is literally labeled `Typical Nucleus (SN)`, and no nucleus
droplet is literally labeled `Typical Cell (SC)`. This literal result is true
by construction: with `assay_col` present, known cells can only be named
`Typical Cell (SC)`, `Nucleus-like Cell (SC)`, or `Neutral`, while known nuclei
can only be named `Typical Nucleus (SN)`, `Cell-like Nucleus (SN)`, or
`Neutral`. Therefore, the useful error counts are the disagreement labels:
59 cells were called nucleus-like (1.97% of all cells), and 7 nuclei were
called cell-like (0.23% of all nuclei).

The signature separates the two known modalities very well when enough
signature evidence is available. On the 5,013 informative droplets, the
continuous `p_cell` score had AUROC 0.9948. At the selected 0.95/0.05 cutoffs,
4,849 droplets received confident calls, of which 4,783 agreed with the known
source: 98.64% accuracy and 98.77% balanced accuracy among confident calls.
Confident-call coverage was 80.82% of all droplets and 96.73% of informative
droplets. The important limitation is unequal coverage: 93.13% for cells but
only 68.50% for nuclei. Of the 1,151 Neutral droplets, 987 had fewer than 10
signature counts, including 873 nuclei, so the lower nucleus coverage mostly
reflects low signature evidence rather than confident cross-modality calls.

Overall conclusion: the independent dataset strongly supports transfer of the
signature—component separation was 0.848 on the sc-fraction scale and the
two-component model was strongly preferred (`BIC1 - BIC2 = 1419.42`). The
signature has excellent discrimination and very few confident errors, but it
abstains frequently on low-count nucleus droplets. Both accuracy and coverage
should be reported; quoting only 98.64% would hide the Neutral calls.

One design caveat is that `suspension_type` and sequencing assay are paired in
this dataset (cells are `10x 5' v1`; nuclei are `10x multiome`). The analysis
therefore demonstrates transfer to this independent mixed dataset, but cannot
separate modality signal from every protocol-specific effect.

### How “how well” was measured

| Quantity | Value | Why report it |
|---|---:|---|
| AUROC on informative droplets | 0.9948 | Threshold-independent separation of known cells and nuclei by `p_cell` |
| Accuracy among confident calls | 98.64% | Fraction of non-Neutral calls agreeing with known source |
| Balanced accuracy among confident calls | 98.77% | Gives cell and nucleus accuracy equal weight |
| Confident-call coverage, all droplets | 80.82% | Shows how often the method makes a call rather than abstaining |
| Confident-call coverage, informative droplets | 96.73% | Separates posterior ambiguity from the low-count gate |
| Cell / nucleus coverage | 93.13% / 68.50% | Exposes the assay-specific coverage imbalance |
| Cross-modality calls | 59 cells; 7 nuclei | Direct answer about confident disagreements |
| Component means | 0.9615 cell; 0.1130 nucleus | Shows biological/statistical separation of the fitted mixture |
| `BIC1 - BIC2` | 1419.42 | Strong support for two components rather than one |

For submission, include the two-row call table above as the primary result and the
diagnostic figure `figures/q1_cellornucem_diagnostics.png`. The table answers
the question exactly; the figure supports it by showing known source versus
prediction on UMAP and the separation of `sc_frac` and `p_cell`. A figure alone
is not sufficient because it does not communicate the Neutral count or exact
error rates.

![Q1a CellorNucEM diagnostics](../figures/q1_cellornucem_diagnostics.png)

### Reproducible run record

- Script: `scripts/q1_mixed.py`
- Command:

  ```bash
  MPLCONFIGDIR=/private/tmp/cellornuc-mpl \
    .venv/bin/python scripts/q1_mixed.py
  ```

- Input: `data/mixed_mBDRC_5pCells_multiomeNuclei_6k.h5ad`
  (`6,000 x 29,987`; 3,000 known cell and 3,000 known nucleus droplets).
- Method: global `Scanpyplus.CellorNucEM` fit using `layer="counts"`,
  `assay_col="suspension_type"`, `sc_label="cell"`,
  `sn_label="nucleus"`, `min_counts=10`, `gamma_hi=0.95`, and
  `gamma_lo=0.05`. These match the mixed-data settings in the supplied demo,
  with label values adapted to this dataset.
- Signature overlap: all 10 sc genes and 9/10 sn genes were present;
  `AC098829.1` was absent.
- Fit diagnostics: cell-component mean sc fraction `0.96145`,
  nucleus-component mean `0.11304`, separation `0.84841`. The two-component
  model was strongly preferred (`BIC1 - BIC2 = 1419.42`).
- Classification counts: 2,735 Typical Cell (SC), 59 Nucleus-like Cell (SC),
  2,048 Typical Nucleus (SN), 7 Cell-like Nucleus (SN), and 1,151 Neutral.
- Independent modality comparison (derived from `p_cell`): 4,849/6,000
  droplets received a confident call (80.82% coverage), and 4,783/4,849 were
  correct (98.64% accuracy among confident calls). Cross-modality calls were
  59/3,000 cells (1.97%) and 7/3,000 nuclei (0.23%). Cell confident-call
  coverage was 93.13%, while nucleus coverage was lower at 68.50%.
- Neutral calls: 987 droplets had fewer than 10 signature counts (114 cells,
  873 nuclei); another 164 informative droplets had posterior probabilities
  between 0.05 and 0.95. This explains most of the lower nucleus coverage.
- Label-semantics caveat: with `assay_col` supplied, a known cell cannot
  literally receive `Typical Nucleus (SN)`, and a known nucleus cannot
  literally receive `Typical Cell (SC)`; those counts are zero by construction.
  The meaningful errors are `Nucleus-like Cell (SC)` and
  `Cell-like Nucleus (SN)`, because they indicate disagreement between the
  unsupervised signature component and the held-out assay identity.
- Outputs: `results/q1/summary.json`, `assay_by_classification.csv`,
  `assay_by_prediction.csv`, `classification_counts.csv`,
  `droplet_classifications.csv.gz`, `mixed_6k_cellornucem.h5ad`, and
  `figures/q1_cellornucem_diagnostics.png`.
