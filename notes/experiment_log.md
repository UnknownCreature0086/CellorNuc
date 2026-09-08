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

### How to interpret each Q1a figure panel

The four panels follow the analysis from known reference labels, to model
calls, to the raw signature statistic, and finally to the fitted posterior.

**Top left — known suspension type on UMAP.** Each point is one droplet. Its
coordinates are the precomputed `X_umap` coordinates supplied in the input
file, and color denotes the known `suspension_type`: cell or nucleus. This is
the reference view against which the prediction panel is compared. The panel
shows that cells and nuclei overlap in several transcriptional neighborhoods,
while some neighborhoods are enriched for one source. It is included to give
biological/cluster context; UMAP is not an input to `CellorNucEM`. UMAP also
does not provide an accuracy estimate because distances are nonlinear,
overplotted points can hide errors, and apparent cluster size/density is not a
confusion matrix.

**Top right — CellorNucEM calls on the same UMAP.** The coordinates and points
are identical to the top-left panel, but colors now show the source-independent
call derived from `p_cell`: `Cell-like`, `Nucleus-like`, or `Neutral`. The EM
fit did not use the known source labels. Broad agreement of the orange and blue
patterns between the two UMAPs is therefore a qualitative indication that the
signature transfers to this dataset. Gray regions reveal where the classifier
abstained, especially in nucleus-rich neighborhoods. This view is useful for
detecting spatially concentrated failures or cell-type-specific patterns, but
the exact 59 cell-to-nucleus-like errors, 7 nucleus-to-cell-like errors, and
1,151 Neutral calls must be read from the call table rather than the UMAP.

**Bottom left — raw cell-signature fraction.** For each informative droplet,
the x-axis is

```text
sc_frac = sc-signature counts / (sc-signature counts + sn-signature counts)
```

and the y-axis is probability density. A value near 1 means almost all measured
signature counts came from the sc gene set; a value near 0 means they came from
the sn gene set. Only droplets with `modality_counts >= 10` are shown so the
ratio is not dominated by extremely sparse evidence. The orange cell
distribution piles up near 1, whereas the blue nucleus distribution is
concentrated near 0 and has a broader right tail. This demonstrates separation
in the actual statistic entering the beta-binomial mixture before posterior
thresholding; the fitted component means, 0.9615 and 0.1130, summarize the same
pattern. Overlap between the curves is where cross-modality or uncertain calls
can arise. This panel does not itself show classification confidence because
the same fraction can carry different evidence when its denominator is 10
versus 1,000.

**Bottom right — fitted posterior probability.** The x-axis is `p_cell`, the
model's posterior probability that a droplet belongs to the high-sc-fraction
component, and the y-axis is density. Color uses known source only for
evaluation after fitting. The cell distribution is concentrated near 1 and the
nucleus distribution near 0, which visually explains the high AUROC of 0.9948.
The dashed vertical lines mark `gamma_lo=0.05` and `gamma_hi=0.95`. For an
informative droplet, values below 0.05 produce a nucleus-like call, values above
0.95 produce a cell-like call, and values in between remain Neutral. Orange
mass on the far left and blue mass on the far right correspond to confident
cross-modality calls. A crucial caveat is that this panel includes low-count
droplets: `CellorNucEM` computes their posterior, but `min_counts=10` prevents
them from receiving a non-Neutral label even if `p_cell` is extreme. Q1b's
count-versus-posterior plot makes that second gate explicit.

Together, the panels support the following reasoning:

```text
known source distribution on UMAP
        -> predicted distribution is qualitatively similar
        -> the raw signature fraction separates the sources
        -> the fitted posterior makes most informative calls highly confident
```

For the final report, a suitable figure caption is: “CellorNucEM generalized
well to independent mixed dataset A. The same UMAP colored by known source and
source-independent model call showed broad concordance, while the raw
cell-signature fraction and fitted `p_cell` posterior separated cells from
nuclei. At 0.05/0.95 posterior cutoffs, 4,783/4,849 confident calls agreed with
source labels (98.64%); 1,151 droplets remained Neutral, primarily because of
the `min_counts` evidence gate analyzed in Q1b.”

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

## 2026-09-07 — Q1b origin of the Neutral calls

### Submission-ready answer

`modality_counts` is the amount of evidence available to `CellorNucEM`, not
the droplet's total UMI count. The source computes

```text
x = sum of raw counts over the sc-signature genes
y = sum of raw counts over the sn-signature genes
modality_counts = m = x + y
```

In this dataset, this is the sum over the 10 detected sc genes and 9 detected
sn genes. With `min_counts=10`, a droplet is informative only when at least 10
raw counts in total occur across those signature genes. A droplet can therefore
have a reasonable library size but still fail this gate if few counts land in
the short signature lists.

The source-code logic explains the Neutral calls precisely. All labels begin
as `Neutral`. A mixed-dataset droplet is overwritten with a typical or
cross-modality label only if both conditions hold:

1. `modality_counts >= min_counts`; and
2. its posterior is confident: `p_cell > gamma_hi` or
   `p_cell < gamma_lo`.

With this run's `gamma_hi=0.95` and `gamma_lo=0.05`, an informative droplet
with `0.05 <= p_cell <= 0.95` remains Neutral. The comparisons in the API are
strict, so equality also remains Neutral. An unrecognized value in
`assay_col` would be a third possible source of Neutral labels, but this
dataset contains only the two recognized values `cell` and `nucleus`.

#### Distribution of `modality_counts`

| Known source | Median | IQR | Mean | Below 10 | Percent below 10 |
|---|---:|---:|---:|---:|---:|
| Cell | 124 | 58–246 | 202.96 | 114/3,000 | 3.8% |
| Nucleus | 24 | 8–59.25 | 50.27 | 873/3,000 | 29.1% |

The nucleus distribution is strongly shifted toward lower signature evidence.
Its first quartile is 8, already below `min_counts=10`, whereas the cell first
quartile is 58. Consequently, nuclei account for 873 of the 987 droplets that
fail the count gate.

#### Exact decomposition of Neutral calls

| Known source | Below `min_counts` | Intermediate posterior | All Neutral |
|---|---:|---:|---:|
| Cell | 114 | 92 | 206 |
| Nucleus | 873 | 72 | 945 |
| **Total** | **987** | **164** | **1,151** |

Thus, 987/1,151 Neutral calls (85.75%) come directly from
`modality_counts < 10`. The remaining 164 (14.25%) pass the count gate but
fall in the intermediate posterior interval. There are no Neutral calls from
unrecognized assay labels and no unexplained states.

The probable technical/biological reasons are:

- `modality_counts` covers only a short signature, so sampling noise and gene
  dropout can leave fewer than 10 relevant counts even when total UMIs are
  higher.
- Nuclei capture less cytoplasmic RNA, while several sc-signature genes are
  abundant cytoplasmic/housekeeping transcripts. The nucleus libraries here
  also come from a different RNA protocol (`10x multiome`) than the cells
  (`10x 5' v1`), which can contribute lower RNA depth.
- One sn-signature gene, `AC098829.1`, is absent from the feature index,
  slightly reducing the measurable evidence available to the model.
- Cell-type, donor, and biological heterogeneity can weaken or mix the two
  signatures. Ambient RNA, damaged cells/nuclei, cytoplasmic leakage, or
  doublets are plausible reasons for intermediate `p_cell`, but this analysis
  does not prove which mechanism applies to an individual droplet.
- The conservative 0.05/0.95 posterior cutoffs intentionally create a broad
  abstention region. Less stringent cutoffs would reduce Neutral calls but
  would trade away confidence; `min_counts` should not be lowered solely to
  force sparse droplets into a class.

Overall, most Neutral calls are deliberate low-evidence abstentions rather
than contradictory classifications. The large excess among nuclei is
explained primarily by their lower `modality_counts` distribution.

### How the distribution was inspected

The figure combines four complementary views:

1. A histogram of `log10(modality_counts + 1)` shows the full distribution,
   including zeros and the long upper tail.
2. An empirical cumulative distribution (ECDF) makes the fraction below the
   threshold directly visible without choosing histogram bins.
3. A stacked count plot gives the exact Neutral-cause decomposition.
4. A `p_cell` versus `modality_counts` plot shows the two independent gates:
   the vertical count boundary and horizontal posterior boundaries.

#### How to interpret each panel in the figure

**Top left — distribution of signature evidence.** The x-axis is
`log10(modality_counts + 1)`, and the y-axis is density rather than the raw
number of droplets. Adding 1 keeps zero-count droplets in the plot, and the log
transformation prevents the long upper tail from compressing all low-count
observations into a narrow area. Density normalization makes the shapes of the
cell and nucleus distributions comparable even if group sizes differ. The
dashed vertical line marks `min_counts=10` on the transformed axis; droplets
to its left fail the evidence gate. The blue nucleus curve is shifted left of
the orange cell curve, showing that nuclei generally contribute fewer counts
from the signature genes. This panel is best for seeing the overall shape,
spread, overlap, and long tail, but exact below-threshold percentages should be
read from the table rather than estimated from histogram area because the
appearance depends on bin placement.

**Top right — empirical cumulative distribution (ECDF).** For any x value,
the y-axis gives the fraction of that source with `modality_counts` at or below
x. The count axis is shown as `modality_counts + 1` on a log scale so zeros can
be included. The dashed line again marks the cutoff. Immediately to the left
of that boundary, the curves correspond to the fractions with fewer than 10
counts: 3.8% of cells and 29.1% of nuclei. The nucleus curve rises earlier and
remains to the left, confirming that low signature evidence is common across
the nucleus distribution rather than being caused by just a few outliers. An
ECDF is useful here because it uses every observation, requires no histogram
bins, and allows a threshold fraction to be read directly.

**Bottom left — exact origin of Neutral calls.** Each bar represents all
Neutral droplets from one known source. Bar height is the number of Neutral
droplets, and colored segments partition that count into mutually exclusive
source-code causes. Gray denotes `modality_counts < 10`; orange denotes
sufficient counts but an intermediate posterior. The nucleus bar is much
taller (945 versus 206), and 873/945 nucleus Neutrals are gray. Across both
sources, the gray segments total 987/1,151, or 85.75% of all Neutrals. This is
the most direct panel for answering *where the Neutral labels came from*. It
shows counts rather than rates, so the accompanying table should be used when
comparing percentages.

**Bottom right — the two classification gates.** Each point is a Neutral
droplet. The x-axis is `modality_counts + 1` on a log scale and the y-axis is
the fitted posterior probability `p_cell`. The vertical dashed line separates
droplets with insufficient signature counts from informative droplets. The
horizontal dotted lines mark `gamma_lo=0.05` and `gamma_hi=0.95`; the shaded
middle band is the posterior abstention interval. Points left of the count
boundary remain Neutral regardless of their posterior, because sparse evidence
is not allowed to produce a class. Points to the right remain Neutral because
their posterior lies between the confidence thresholds. This panel explains
the logical AND in the source code: a non-Neutral label requires both enough
counts and a posterior outside the middle band. The plotted posterior for a
low-count droplet should not be interpreted as a reliable call—the API
computes it, but deliberately refuses to classify from it.

Together, the panels move from descriptive evidence to mechanism:

```text
different count distributions
        -> many nuclei fall below min_counts
        -> low-count droplets dominate the Neutral category
        -> remaining Neutrals are posterior-based abstentions
```

For the final report, the most defensible one-sentence figure interpretation
is: “Nuclei had a markedly left-shifted signature-count distribution, causing
29.1% to fail the `min_counts=10` evidence gate; accordingly, low signature
counts explained 987/1,151 (85.75%) Neutral calls, while only 164 informative
droplets were Neutral because their posterior lay between 0.05 and 0.95.”

For a submitted answer, include the distribution/threshold table, the Neutral
decomposition table, and the figure below. The ECDF or histogram establishes
the distributional shift; the exact tables prevent the visual from being
overinterpreted.

![Q1b origins of Neutral calls](../figures/q1b_neutral_origins.png)

### Reproducible run record

- Script: `scripts/q1b_neutral.py`
- Command:

  ```bash
  MPLCONFIGDIR=/private/tmp/cellornuc-mpl \
    .venv/bin/python scripts/q1b_neutral.py
  ```

- The script reruns the same global `CellorNucEM` API call used for Q1a and
  mirrors the public API's count and posterior gates when assigning an
  explanatory `neutral_cause`.
- Outputs:
  - `results/q1/q1b_modality_counts_summary.csv`
  - `results/q1/q1b_modality_counts_bins.csv`
  - `results/q1/q1b_neutral_causes.csv`
  - `results/q1/q1b_neutral_droplets.csv.gz`
  - `results/q1/q1b_summary.json`
  - `figures/q1b_neutral_origins.png`

## 2026-09-07 — Q2 single-modality sc dataset B (6k)

### Submission-ready answer

I ran a global `CellorNucEM` fit on dataset B's raw `counts` layer using all
10 supplied sc-signature genes and all 10 supplied sn-signature genes. Because
the file is known to contain an sc library but has no assay column, I added the
constant metadata column `modality="sc"`. As specified in the source, this
column is used only to name the final classes; it is not used to fit the
mixture. I followed the supplied single-modality demo with `min_counts=10`,
`gamma_lo=0.05`, `gamma_hi=0.95`, `max_p_nuc=0.2`, and `min_p_cell=0.8`.

**Q2a.** CellorNucEM classified **180 of 6,000 droplets as `Nucleus-like Cell
(SC)`**, a fraction of **0.0300, or 3.00%**.

| CellorNucEM classification | Droplets | Fraction | Percent |
|---|---:|---:|---:|
| Typical Cell (SC) | 5,820 | 0.9700 | 97.00% |
| Nucleus-like Cell (SC) | 180 | 0.0300 | 3.00% |

**Q2b.** The `sc_frac` histogram is **bimodal**, but highly imbalanced. It has
a dominant cell-like mode near 1 and a much smaller nucleus-like mode near 0,
with a broad bridge of observations between them. In the 40-bin histogram,
2,583 observations fall in `[0.975, 1.000]`, while 265 fall in `[0, 0.025)`.
The edge spikes should not be read as perfectly measured latent fractions:
exact ratios of 0 or 1 are especially common when only a modest number of
signature counts were observed.

![Q2 sc_frac histogram](../figures/q2_sc_frac_histogram.png)

### What the numbers mean

For each droplet, the implementation computes

```text
x = raw counts over the sc-enriched genes
y = raw counts over the sn-enriched genes
m = modality_counts = x + y
sc_frac = x / m
```

Thus, `sc_frac` near 1 means that nearly all observed signature counts came
from the sc-enriched genes; `sc_frac` near 0 means that nearly all came from
the sn-enriched genes. It is a signature-composition statistic, not the
fraction of the whole transcriptome that is cytoplasmic and not itself a
posterior probability.

Of the 6,000 droplets, 5,973 have defined `sc_frac`; the remaining 27 have
`m=0`, so their ratio is undefined. The overall median defined value is
0.9579. The 180 nucleus-like hard calls have median `sc_frac=0`, mean 0.0247,
and range 0–0.125. By contrast, the droplets ultimately labeled typical have
median 0.9625. The histogram therefore exposes a distinct low-fraction tail
that a single summary such as the median would hide.

The constrained beta-binomial mixture estimated component mean fractions of
0.9091 for the cell-like component and 0.2000 for the nucleus-like component,
giving a separation of 0.7091. The two-component BIC was lower than the
one-component BIC by 662.51, which is strong statistical support for two
distributional components in the informative observations. However, the
nucleus mean equals its imposed upper bound of 0.2. It must therefore be read
as a constrained anchor, not as a freely estimated biological constant.
Neither the BIC nor the two-component fit proves that every low-mode droplet
is a lysed cell; cell-type-specific signature behavior and technical effects
remain alternatives.

### Continuous score versus hard classification

The source converts `(x,m)` into `p_cell`, the fitted posterior probability of
membership in the high-`sc_frac` component. It then calls an sc droplet
`Nucleus-like Cell (SC)` only when both conditions are true:

```text
modality_counts >= 10  AND  p_cell < 0.05
```

The decision accounting is:

| Model/evidence state | Droplets | Percent of all droplets | Final sc-only label |
|---|---:|---:|---|
| Confident nucleus-like (`p_cell < 0.05`) | 180 | 3.00% | Nucleus-like Cell (SC) |
| Intermediate posterior (`0.05 <= p_cell <= 0.95`) | 508 | 8.47% | Typical Cell (SC) |
| Confident cell-like (`p_cell > 0.95`) | 4,600 | 76.67% | Typical Cell (SC) |
| Below `min_counts` | 712 | 11.87% | Typical Cell (SC) |

This table explains a subtle source-code behavior: single-modality data have
no displayed Neutral class. All 1,220 droplets that would be undecided because
of intermediate posterior or low evidence are folded back into their known
assay's typical label. Consequently, `Typical Cell (SC)` does **not** mean
that every one of its 5,820 members had `p_cell > 0.95`.

The raw ratio and posterior are connected but not interchangeable. `p_cell`
uses the fitted component shapes, their mixture weights, and the denominator
`m`. Two droplets with the same `sc_frac` can have different posterior
certainty when one has 10 signature counts and another has 1,000. The
posterior thresholds then turn that continuous evidence into conservative
hard labels, while `min_counts` prevents a sparse ratio from becoming a hard
call. The second figure makes these layers visible.

![Q2 score to classification](../figures/q2_score_to_classification.png)

The fitted cell mixture weight was 0.93525, corresponding to a soft
nucleus-component weight of 0.06475 among informative droplets. Summing each
informative droplet's posterior nucleus probability gives 342.06 expected
nucleus-component memberships, or 6.47% of the 5,288 informative droplets.
This is deliberately larger than the 3.00% hard-call fraction: soft membership
credits ambiguous droplets fractionally, whereas the hard rate counts only
droplets past the stringent 0.05 cutoff and the count gate. For the quiz's
literal Q2a wording, 3.00% is the answer; 6.47% is a useful model diagnostic,
not a replacement for it.

### Why the fraction matters, and its limitations

The 3% fraction is an estimate of the library's confidently flagged
nucleus-like burden. In a real sc project it can guide QC, flag libraries or
cell populations for review, and quantify how much downstream expression,
clustering, differential-expression, or abundance analysis might be affected
by droplets with depleted cytoplasmic signal. A higher fraction would raise
greater concern about preparation quality or contamination. The scientific
cost of removing 3% may still be nontrivial if those droplets are concentrated
in a rare cell type.

They are not distributed uniformly across the supplied annotations. For
example, 35/69 annotated `CCD-PC` droplets and 15/25 `CCD-PC|CNT doub`
droplets were flagged, while many larger cell-type groups had much lower
rates. These small-group percentages are imprecise, but the concentration is
a warning against interpreting every flag as proven cytoplasmic lysis or
blindly filtering all 180 droplets. It could reflect true cell-type-specific
susceptibility, annotation/doublet effects, or lack of universality of the
short signature. The saved cell-type and batch tables support that review;
orthogonal QC evidence would be needed to establish mechanism.

In general:

- a clearly **unimodal high-`sc_frac`** distribution would suggest one largely
  intact-cell population and little evidence for a distinct nuclear-like
  subgroup;
- a clearly **bimodal** distribution suggests two compositional regimes and
  makes a cell-like/nucleus-like mixture interpretation plausible;
- modes alone are not classifications. Their overlap, count depth, fitted
  posterior, and thresholds determine which individual droplets receive hard
  calls.

### Reproducible run record

- Script: `scripts/q2_single.py`
- Command:

  ```bash
  MPLCONFIGDIR=/private/tmp/cellornuc-mpl \
    XDG_CACHE_HOME=/private/tmp/cellornuc-cache \
    .venv/bin/python scripts/q2_single.py
  ```

- Input: `data/single_KidneyRaji_sc_6k.h5ad` (`6,000 x 28,562`).
- Method: global constrained `CellorNucEM` on `layer="counts"`; a constant
  `modality="sc"` was added for post-fit class naming.
- Signature overlap: 10/10 sc genes and 10/10 sn genes.
- Outputs:
  - `results/q2/summary.json`
  - `results/q2/classification_counts.csv`
  - `results/q2/decision_gate_counts.csv`
  - `results/q2/sc_frac_histogram_bins.csv`
  - `results/q2/sc_frac_summary.csv`
  - `results/q2/nucleus_like_by_cell_type.csv`
  - `results/q2/nucleus_like_by_batch.csv`
  - `results/q2/droplet_classifications.csv.gz`
  - `results/q2/single_6k_cellornucem.h5ad`
  - `figures/q2_sc_frac_histogram.png`
  - `figures/q2_score_to_classification.png`
