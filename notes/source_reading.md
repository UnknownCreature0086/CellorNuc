# Source reading: `CellorNuc` and `CellorNucEM`

This note follows `scanpyplus/Scanpyplus.py` lines 1351-1841. The quiz uses
`CellorNucEM`, but the older `CellorNuc` function is useful as a contrast.

## 1. Mental model

The EM classifier reduces every droplet to counts from two short gene sets:

```text
AnnData raw counts
    |
    +-- x = sum of counts for sc-enriched genes
    +-- y = sum of counts for sn-enriched genes
    |
    +-- m = x + y                 evidence carried by the signatures
    +-- sc_frac = x / m           observed cell-like fraction
    |
    +-- fit two beta-binomial components to (x, m)
            |
            +-- high-mean component = cell-like
            +-- low-mean component  = nucleus-like
            |
            +-- p_cell = posterior probability of high-mean component
                    |
                    +-- threshold p_cell
                    +-- combine with optional assay metadata to name the class
```

The method is therefore a classifier of **cell-like versus nucleus-like
signature composition**. It does not use all genes in the likelihood, and
`modality_counts` is not the droplet's total UMI count: it is the sum across
the sc and sn signature genes only.

## 2. Vocabulary and notation

| Symbol/code name | Meaning |
|---|---|
| `sc_genes` | Genes whose expression is enriched in intact-cell data. |
| `sn_genes` | Genes whose expression is enriched in nucleus data. |
| `x` | Per-droplet sum over the sc signature genes. |
| `y` | Per-droplet sum over the sn signature genes. |
| `m = x + y` | Total signature evidence for that droplet. |
| `sc_frac = x / m` | Observed fraction of signature counts coming from sc genes. |
| `(a, b)` or `ab` | Shape parameters of a fitted beta distribution. |
| `p = a / (a + b)` | Mean sc fraction of one component. |
| `lam_cell` | Estimated fraction of observations in the cell-like component. |
| `gamma` / `p_cell` | Posterior probability that a droplet belongs to the cell-like component. |
| `enough` | Whether `m >= min_counts`; only these droplets fit or receive evidence-based labels. |

`sc` and `sn` describe the known assay/source modality. “Cell-like” and
“nucleus-like” describe what the signature model sees. That is why a label such
as `Nucleus-like Cell (SC)` is meaningful: the metadata says the droplet came
from an sc assay, while the fitted signature component looks nuclear.

## 3. Gene signatures and `LoadGeneSignatures()` (lines 1351-1391)

### `gene_signatures.json`

The JSON file currently has one signature set, `hcka_v1`, containing ten
sc-enriched and ten sn-enriched genes. It also stores provenance rather than
only gene names:

- the signatures were derived from Human Kidney Cell Atlas sc and sn objects;
- mitochondrial and ribosomal genes were excluded during feature selection;
- several feature-ranking/statistical signals were cross-referenced;
- the set is marked `provisional - pre statistical validation`.

That last point is important scientifically: the JSON records a useful trained
signature, not a universal biological law. The quiz's independent-data check
is therefore a validation exercise.

Lines 1351-1352 also define module-level `sc_genes` and `sn_genes` lists that
duplicate the current JSON values. `LoadGeneSignatures()` does not read those
variables; it reads the JSON file.

### Loader, block by block

- **1354-1378:** The signature name defaults to `hcka_v1`. Passing
  `signature=None` requests the entire metadata dictionary instead of a pair of
  gene lists. `path` can point to a different signature file.
- **1379-1382:** If no path is supplied, it builds a path relative to
  `Scanpyplus.py` itself. This makes loading independent of the current working
  directory.
- **1383-1384:** Standard `json.load` converts the file into nested Python
  dictionaries/lists.
- **1385-1386:** `signature=None` returns everything, including provenance.
- **1387-1389:** An unknown signature name raises a helpful error listing the
  available keys.
- **1390-1391:** A named signature returns exactly two lists:
  `(sig['sc_genes'], sig['sn_genes'])`.

The ordinary entry point is:

```python
sc_genes, sn_genes = Scanpyplus.LoadGeneSignatures()
```

## 4. Older score-based `CellorNuc` (lines 1393-1503)

This function and `CellorNucEM` are not two names for the same algorithm.

### Inputs and data assumptions (1393-1441)

`CellorNuc` requires:

1. an `AnnData` object;
2. both signature lists;
3. an `obs` cell-type column;
4. an `obs` assay/modality column with recognizable sc and sn labels.

Unlike the EM implementation, it calls `scanpy.tl.score_genes` on the active
matrix. That score is normally interpreted on appropriately normalized
expression data. The function does not select a layer or verify the matrix's
scale.

`copy=False` means mutate the supplied object and return `None`.
`copy=True` means work on and return an independent `AnnData` object.

### Execution (1443-1503)

- **1443:** Apply the copy policy.
- **1445-1449:** Require both requested `obs` columns.
- **1451-1456:** Intersect both signature lists with `adata.var_names`; fail if
  either intersection is empty.
- **1458-1460:** Use Scanpy's gene scoring to create `sc_score` and `sn_score`.
- **1462-1466:** Define `cell_vs_nucleus_score = sc_score - sn_score`, then
  create a percentile-clipped version.
- **1468-1469:** Initialize every label as `Uncertain`.
- **1471-1487:** Within each cell type, use known sc observations to estimate a
  lower sc boundary and known sn observations to estimate an upper sn boundary.
  The boundaries are forced to opposite sides of zero.
- **1489-1497:** Compare each known sc/sn droplet with those boundaries to call
  typical or cross-modality profiles.
- **1499-1503:** Convert labels to a categorical column and obey the copy/return
  policy.

This method is metadata-calibrated: known assay labels participate in estimating
its thresholds. By contrast, `CellorNucEM` fits signature components without
using `assay_col`; that column only changes the final names.

### Two source-level caveats

1. `cell_vs_nucleus_score_clipped` is stored at lines 1465-1466, but lines
   1479-1497 classify with the **unclipped** `cell_vs_nucleus_score`. Thus
   `clip_percentiles` does not affect the calls in this implementation.
2. Lines 1491-1497 assign strings such as `Typical Cell (SC)`, but line 1500's
   categorical levels omit the `(SC)`/`(SN)` suffixes. Pandas converts assigned
   strings that are absent from the declared categories to missing values. This
   should be treated as an implementation inconsistency, not intended biology.

For the starter quiz, the newer EM function is the relevant implementation.

## 5. `_set_counts_EM` (lines 1508-1514)

This helper converts a signature list into one count per droplet.

- **1509:** Keep only genes present in `adata.var_names`. Partial overlap is
  allowed.
- **1510-1511:** Zero overlap is fatal because the corresponding signal cannot
  be measured.
- **1512:** Convert gene names to integer column positions.
- **1513:** Use `adata.layers[layer]` when a layer is named; otherwise use
  `adata.X`.
- **1514:** Subset to the signature columns, sum across genes for each droplet,
  and flatten the result into a one-dimensional NumPy array.

Calling it once on `sc_genes` produces `x`; calling it on `sn_genes` produces
`y`. Because the beta-binomial likelihood is a count model, this is why the
demo emphasizes `layer='counts'`.

## 6. Beta-binomial mathematics (lines 1517-1556)

### What counts as a “trial” and a “success” here?

For droplet $i$:

- $x_i$ is its count across sc-enriched signature genes;
- $y_i$ is its count across sn-enriched signature genes;
- $m_i=x_i+y_i$ is its total signature count;
- an sc-signature count is treated as a “success,” and an sn-signature count as
  a “failure.”

The likelihood is conditional on $m_i$. It models how the fixed amount of
signature evidence divides between the two gene sets; it does not model the
droplet's total transcript count or why $m_i$ varies between droplets.

This is a working statistical abstraction. Real UMIs from different genes are
not literally identical independent coin tosses. The beta-binomial summarizes
the observed composition and its extra heterogeneity without trying to model
every gene separately.

### The hierarchical model

A binomial model would assign one fixed sc fraction $p$ to every droplet in a
component. Biological and technical heterogeneity make that too restrictive.
The beta-binomial instead gives each droplet a latent fraction $q_i$:

```text
q_i | component k  ~ Beta(a_k, b_k)
x_i | q_i, m_i, k  ~ Binomial(m_i, q_i)
```

Integrating out $q_i$ gives

\[
X_i\mid m_i,k \sim \operatorname{BetaBinomial}(m_i,a_k,b_k).
\]

Its probability mass function is

\[
f_k(x_i\mid m_i)
= {m_i \choose x_i}
  \frac{B(x_i+a_k,\;m_i-x_i+b_k)}{B(a_k,b_k)},
\]

where $B(\cdot,\cdot)$ is the beta function. This is what
`scipy.stats.betabinom.logpmf(x, m, a, b)` evaluates in log form.

### Statistical meaning of `a` and `b`

The raw shape parameters are easiest to interpret through two derived values:

\[
p_k=\frac{a_k}{a_k+b_k}, \qquad
\kappa_k=a_k+b_k.
\]

- $p_k$, the component mean, is the expected sc-signature fraction.
- $\kappa_k$, the concentration, controls how tightly droplets cluster around
  that mean.

The latent fraction has

\[
E(q_i\mid k)=p_k,
\qquad
\operatorname{Var}(q_i\mid k)
=\frac{p_k(1-p_k)}{\kappa_k+1}.
\]

Thus two components can have the same mean but different heterogeneity. A large
$\kappa$ creates a narrow Beta distribution and approaches an ordinary
binomial; a small $\kappa$ permits substantial droplet-to-droplet variation.
Neither `a` nor `b` is a number of cells, genes, or observed UMIs. They are
fitted distribution-shape parameters.

For observed counts,

\[
E(X_i\mid m_i,k)=m_i p_k,
\]

\[
\operatorname{Var}(X_i\mid m_i,k)
=m_i p_k(1-p_k)\frac{m_i+\kappa_k}{1+\kappa_k}
=m_i p_k(1-p_k)\left[1+(m_i-1)\rho_k\right],
\]

with intra-class correlation $\rho_k=1/(\kappa_k+1)$. The bracketed factor is
the extra variance relative to a binomial. This overdispersion is the main
reason to use beta-binomial rather than binomial.

For the observed ratio $X_i/m_i$ when $m_i>0$,

\[
\operatorname{Var}(X_i/m_i\mid m_i,k)
=\frac{p_k(1-p_k)(m_i+\kappa_k)}{m_i(1+\kappa_k)}.
\]

Small $m_i$ therefore makes `sc_frac` noisy. Increasing $m_i$ reduces
binomial sampling noise, but the variance does not collapse to zero because the
Beta layer represents persistent between-droplet heterogeneity.

### Code-to-math dictionary

| Code value | Mathematical role | Shape/interpretation |
|---|---|---|
| `x`, `y`, `m` | $x_i,y_i,m_i$ | Length `n_obs` arrays, one value per droplet. |
| `enough` | $1[m_i\geq\text{min_counts}]$ | Boolean gate for fitting and confident classification. |
| `pairs[:, 0]`, `xu` | Unique $x_j$ values | sc counts for each unique `(x,m)` pattern. |
| `pairs[:, 1]`, `mu` | Unique $m_j$ values | Total signature counts for each unique pattern; `mu` does not mean a statistical mean here. |
| `cnt`, `w` | $c_j$ or soft weight | Pattern frequency, possibly multiplied by an EM responsibility. |
| `p` in line 1583 | $x_j/m_j$ | Raw observed fraction used only to initialize groups. |
| `ab1`, `ab_cell` | $(a_C,b_C)$ | Shape parameters for the high-mean/cell-like component. |
| `ab0`, `ab_nuc` | $(a_N,b_N)$ | Shape parameters for the low-mean/nucleus-like component. |
| `lam`, `lam_cell` | $\lambda=P(Z=C)$ | Mixture prevalence/prior probability of the cell-like component. |
| `l1`, `l0` | log joint component scores | Log mixture weight plus log beta-binomial PMF. |
| `norm` | log marginal probability | `logsumexp(l1,l0)` for each unique pattern. |
| `r1` | $r_j=P(Z_j=C\mid x_j,m_j)$ | E-step cell responsibility for each unique pattern. |
| `ll`, `prev_ll` | $\ell$ | Frequency-weighted mixture log likelihood. |
| `gamma`, `p_cell` | $P(Z_i=C\mid x_i,m_i)$ | Final posterior for every droplet. |

There are three different quantities that may all be informally called “p”:

1. `p = xu/mu` is an observed ratio used during initialization;
2. $p_C=a_C/(a_C+b_C)$ and $p_N=a_N/(a_N+b_N)$ are component means;
3. `p_cell` is the posterior probability of component membership.

Keeping these distinct is essential when reading the code or answering a
question about “component means p.”

### Weighted maximum likelihood: `_fit_bb_weighted_EM`

For one component, the weighted log likelihood is

\[
\ell(a,b)=\sum_j w_j\log f(x_j\mid m_j,a,b).
\]

When `w=cnt`, this is exactly the ordinary likelihood of all droplets after
identical `(x,m)` observations have been compressed. For example, a unique pair
seen 30 times contributes 30 copies of the same log PMF. During EM, `w` becomes
`cnt*r1` or `cnt*(1-r1)`, producing fractional/soft component membership.

There is no simple closed-form MLE for beta-binomial `a,b`, so lines 1530-1542
minimize the negative weighted log likelihood numerically:

- `log_ab` is optimized and exponentiated, guaranteeing $a,b>0$;
- `init=(5,5)` is only the optimizer's starting point, not a Bayesian prior;
- invalid PMFs and concentration above `max_conc=1e8` receive a large penalty;
- Nelder-Mead returns fitted shape parameters.

The concentration cap prevents numerical runaway when a component is almost
binomial. As $\kappa\to\infty$, between-droplet heterogeneity vanishes, so very
large concentrations are practically indistinguishable while becoming
numerically unstable.

The optimizer result's success flag is not checked. Fitted means, separation,
BIC, convergence messages, and sensitivity to initialization should therefore
be treated as diagnostics rather than assuming every returned tuple is ideal.

### `_mean_p_EM` and `_constrain_EM`

`_mean_p_EM((a,b))` computes $p=a/(a+b)$.

For concentration $\kappa=a+b$, `_constrain_EM` clips the mean to a target
$p^*$ and returns

\[
(a^*,b^*)=(p^*\kappa,(1-p^*)\kappa).
\]

This changes the component location while preserving its concentration. It is
a post-M-step projection, not a fresh numerical maximization of the likelihood
subject to the bound. Consequently, with active constraints the procedure is
best understood as projected/constrained EM, and the usual guarantee that an
unconstrained EM likelihood never decreases need not apply to the projection.

`anchor_lo`/`anchor_hi` do not call this function and are not constraints.

## 7. Two-component fit: `_fit_two_comp_EM` (lines 1559-1637)

The mixture has a high-fraction component intended to represent cell-like
profiles and a low-fraction component intended to represent nucleus-like
profiles.

Introduce an unobserved indicator $Z_i\in\{C,N\}$:

\[
P(Z_i=C)=\lambda,\qquad P(Z_i=N)=1-\lambda.
\]

Conditional on that indicator, each component has its own beta-binomial shape:

\[
X_i\mid m_i,Z_i=C\sim BB(m_i,a_C,b_C),
\]

\[
X_i\mid m_i,Z_i=N\sim BB(m_i,a_N,b_N).
\]

Because $Z_i$ is not observed by the fitting algorithm, the marginal PMF is

\[
P(X_i=x_i\mid m_i)
=\lambda f_C(x_i\mid m_i)+(1-\lambda)f_N(x_i\mid m_i).
\]

For compressed patterns $j$ with frequencies $c_j$, the quantity optimized
by the code is

\[
\ell(\theta)=
\sum_j c_j\log\left[
\lambda f_C(x_j\mid m_j)+(1-\lambda)f_N(x_j\mid m_j)
\right],
\]

where
$\theta=(a_C,b_C,a_N,b_N,\lambda)$. This is a five-parameter mixture.

### Arguments and return value (1559-1581)

- `xu`, `mu`, `cnt`: unique `(x, m)` pairs and their frequencies.
- `init_from`: an existing fit used as a warm start, especially the global fit
  before fitting individual cell types.
- `init`: initialization scheme; documented choices are `anchored` and `split`.
- `anchor_lo`, `anchor_hi`: tails used only by anchored initialization.
- `max_p_nuc`, `min_p_cell`: optional hard bounds applied after EM updates.
- `max_iter`, `tol`: stopping controls.

It returns the two `(a, b)` pairs, `lam_cell`, mixture log likelihood, and BIC
values for one- and two-component models.

### Initialization precedence (1582-1603)

1. **Warm start, 1584-1586:** If `init_from` is supplied, its component
   parameters and mixing weight take precedence over `init`.
2. **Anchored, 1587-1595:** Compute `p=x/m`; initialize the cell component from
   unique pairs above `anchor_hi` and the nucleus component from pairs below
   `anchor_lo`. If a tail has fewer than two unique `(x,m)` pairs, use weak
   fallback shapes `(19,1)` (mean 0.95) or `(1,19)` (mean 0.05). Initialize the
   cell mixture weight from the high-tail observation count, clipped to
   `[0.01, 0.99]`; if there is no high-tail pair at all, it defaults to `0.9`.
3. **Split, 1596-1603:** Split at the evidence-weighted average fraction. Since
   the weights are `cnt * m`, the split equals the aggregate sc-signature count
   divided by aggregate signature count. Fit one component on each side. If
   either side has fewer than two unique pairs, initialization fails.

The source uses a plain `else`, so any `init` value other than the exact string
`"anchored"` follows the split branch. A misspelling is therefore not rejected.

The fallback tuples `(19,1)` and `(1,19)` are sometimes described in the code
comment as “weak priors,” but mathematically they are **initial parameter
values**, not priors used in Bayesian inference. Once EM starts, the optimizer
can move them unless the separate mean constraints are active.

For split initialization, the code's threshold has a useful simplification:

\[
\text{split}
=\frac{\sum_j(c_jm_j)(x_j/m_j)}{\sum_j c_jm_j}
=\frac{\sum_j c_jx_j}{\sum_j c_jm_j}.
\]

It is the aggregate fraction of all signature counts assigned to sc genes, not
the unweighted average droplet fraction.

### `anchored` versus `split`

| Property | `init="anchored"` | `init="split"` |
|---|---|---|
| Initial groups | Fixed low and high tails of `sc_frac` | Two sides of the data's weighted mean |
| Default cutoffs | `< 0.3` and `> 0.7` | Dataset-dependent |
| Empty-tail behavior | Uses a weak 0.05/0.95 fallback component | Raises if a side has too few unique pairs |
| Intended advantage | Keeps initial components at biologically distinct ends under imbalance | Lets both starts be determined by the observed distribution |
| Main risk | Tail thresholds may not match a shifted dataset | In single-modality data, the mean may cut the dominant peak into two pieces |

Both are only starting strategies. Ordinary EM can move the component means
away from their starts. To keep them within fixed biological ranges throughout
fitting, use `max_p_nuc` and/or `min_p_cell`.

Mixture likelihoods are generally non-convex and can have local optima. Two
initializations can therefore converge to different parameter sets even though
they use exactly the same observations and likelihood. This sensitivity is
most relevant when components overlap strongly, one component is rare, or the
data are effectively single-modality.

### EM iterations (1605-1625)

Each loop alternates between two steps:

**E-step (1607-1612)**

Given current parameters, the responsibility of the cell-like component is

\[
r_j=P(Z_j=C\mid x_j,m_j)
=\frac{\lambda f_C(x_j\mid m_j)}
       {\lambda f_C(x_j\mid m_j)+(1-\lambda)f_N(x_j\mid m_j)}.
\]

The code first constructs

\[
\ell_{Cj}=\log\lambda+\log f_C(x_j\mid m_j)
\]

and

\[
\ell_{Nj}=\log(1-\lambda)+\log f_N(x_j\mid m_j)
\]

as `l1` and `l0`. `logsumexp([l1,l0])` is the stable log denominator. Then
`r1 = exp(l1 - norm)`. `ll = sum(cnt*norm)` is the observed-data mixture log
likelihood, not the complete-data likelihood.

`r1` is a temporary posterior for each **unique pattern** under the current
iteration's parameters. After convergence, `_posterior_cell_EM` applies the
same formula to every original droplet and stores it as `gamma`/`p_cell`.
Droplets with identical `(x,m)` and the same fit necessarily receive identical
posteriors.

**M-step (1615-1625)**

The mixture weight has the closed-form update

\[
\lambda^{new}=\frac{\sum_j c_jr_j}{\sum_jc_j}.
\]

This is why `lam_cell` can be interpreted as the fitted cell-component
prevalence: it is the mean soft assignment, not necessarily the fraction of
droplets that pass the final `gamma_hi` threshold.

The shape updates maximize the expected component-specific log likelihoods:

\[
(a_C,b_C)^{new}
=\arg\max_{a,b}\sum_j c_jr_j\log f(x_j\mid m_j,a,b),
\]

\[
(a_N,b_N)^{new}
=\arg\max_{a,b}\sum_j c_j(1-r_j)\log f(x_j\mid m_j,a,b).
\]

These have no closed-form solution, so each M-step calls the numerical weighted
beta-binomial fitter. The previous shapes are passed as its next starting point.
Afterward, optional mean bounds project the updated shapes as described above.

The loop stops only after at least four iterations when the log-likelihood
improvement is smaller than `tol*abs(ll)`. A negative improvement also satisfies
that test. Moreover, if the condition fires, the returned `prev_ll` is the
previous accepted value rather than the just-computed `ll`; normally the
difference is tiny, but it is worth knowing when auditing convergence.

### Component identity and BIC (1627-1637)

- **1627-1628:** If necessary, swap components so `ab_cell` is always the
  higher-mean one. “Cell” is therefore defined operationally as the component
  with more sc-signature composition.
- **1630-1632:** Fit a single beta-binomial component to the same data.
- **1633-1637:** Store both models' BIC values. The one-component model counts
  two parameters (`a,b`); the mixture counts five (`a,b` twice plus `lam`).
  Lower BIC is preferred. `bic_2 < bic_1` is evidence that two components fit
  better after penalizing complexity, not proof that they are two biological
  modalities.

The exact formulas in the return dictionary are

\[
BIC_1=-2\ell_1+2\log n,
\qquad
BIC_2=-2\ell_2+5\log n,
\]

where $n=\sum_jc_j$ is the number of informative droplets, not the number of
unique pairs and not their total signature counts. A useful comparison is
$\Delta BIC=BIC_1-BIC_2$: positive values favor the two-component fit, while
negative values favor one component. The code uses only the sign and emits a
warning when $BIC_2\geq BIC_1$.

Mixture models have nonstandard boundary behavior, especially when a component
weight approaches zero or two components coincide. BIC is best used here as a
diagnostic alongside component separation, stability across initializations,
and biological plausibility.

### Identifiability in single-modality data

A two-component model can mathematically divide one broad or skewed population
into two fitted components even when there are not two biological modalities.
In that setting:

- `split` may divide the main peak at its weighted center;
- `anchored` supplies starts at opposite tails but cannot create evidence that
  is absent from the likelihood;
- `max_p_nuc` and `min_p_cell` encode external structural assumptions;
- BIC, fitted separation, component prevalence, and initialization sensitivity
  help reveal whether the minority component is well supported.

The high-mean/low-mean names resolve statistical label switching, but they do
not by themselves establish the biological identity of two components.

## 8. `_posterior_cell_EM` (lines 1640-1645)

This helper applies Bayes' rule to each `(x,m)` observation:

```text
p_cell = P(cell component | x, m)
       = lam_cell * BB(x | m, ab_cell)
         ---------------------------------------------------------
         lam_cell * BB(x | m, ab_cell)
         + (1-lam_cell) * BB(x | m, ab_nuc)
```

It performs the calculation in log space for numerical stability.

The same relationship is especially interpretable on the log-odds scale:

\[
\log\frac{p_{cell,i}}{1-p_{cell,i}}
=\log\frac{\lambda}{1-\lambda}
+\log\frac{f_C(x_i\mid m_i)}{f_N(x_i\mid m_i)}.
\]

The posterior log odds equal **prior/prevalence log odds** plus a **data log
likelihood ratio**. Consequently:

- `lam_cell` shifts every posterior toward the more prevalent fitted component;
- the relative component PMFs measure how compatible the droplet is with each
  component;
- `m` affects how much evidence the observed ratio carries.

This explains why `p_cell` is not interchangeable with `sc_frac`:

- `sc_frac` is just the observed ratio `x/m`;
- `p_cell` also accounts for the fitted component locations and dispersions,
  the number of signature counts `m`, and the fitted component prevalence
  `lam_cell`.

Two droplets with the same `sc_frac` but different `m` can therefore have
different certainty.

For a purely hypothetical illustration, suppose the fitted cell component is
`Beta(8,2)`, the nucleus component is `Beta(2,8)`, and `lam_cell=0.7`. The
following observations all have `sc_frac=0.8`, but their posteriors differ:

| `x` | `m` | `sc_frac` | Approximate `p_cell` |
|---:|---:|---:|---:|
| 4 | 5 | 0.8 | 0.9747 |
| 8 | 10 | 0.8 | 0.9929 |
| 40 | 50 | 0.8 | 0.9995 |

The ratio points in the same direction each time, while more observations make
the likelihood ratio more decisive for these particular component parameters.
This example is explanatory only and is not calculated from either quiz file.

At the extreme $m=0,x=0$, both beta-binomial PMFs equal 1, so the posterior is
just `lam_cell`: there is no droplet-specific signature evidence. The
`enough` gate prevents such a posterior from becoming an evidence-based label.
For $0<m<\text{min_counts}$, the model can still compute a posterior, but the
classification gate treats the evidence as insufficient.

`p_cell=0.9` should be read as “90% posterior membership in the fitted
high-sc-fraction component under this model and its estimated parameters.” It
is not automatically a 90% externally calibrated probability that the droplet
is biologically an intact cell, nor is it a frequentist confidence interval.
That interpretation requires the component-to-biology mapping and model to be
credible on the dataset at hand.

## 9. Main API: `CellorNucEM` (lines 1648-1841)

### Signature and parameters (1648-1675)

| Parameter | Role |
|---|---|
| `adata` | Input `AnnData`; observations are droplets and variables are genes. |
| `sc_genes`, `sn_genes` | Signature lists, usually returned by `LoadGeneSignatures`. |
| `assay_col` | Optional metadata used **only after fitting** to name typical/anomalous classes. |
| `cell_type_col` | Optional metadata for separate within-cell-type refits. |
| `sc_label`, `sn_label` | Exact values in `assay_col` that mean sc and sn. |
| `layer` | Raw-count layer; `None` uses `adata.X`. |
| `min_counts=10` | Minimum signature count `m`, not minimum whole-transcriptome UMI count. |
| `gamma_hi`, `gamma_lo` | Posterior thresholds for confident cell-like and nucleus-like calls. |
| `min_cells_fit` | Minimum informative observations needed before trying a cell-type-specific fit. |
| `min_separation` | Minimum difference between fitted component means for accepting a cell-type fit. It is not applied to the global fit. |
| `init`, `anchor_lo`, `anchor_hi` | Global initialization behavior. |
| `max_p_nuc`, `min_p_cell` | Optional constraints on fitted component means. |
| `verbose` | Print every EM iteration. |
| `print_mode` | Print summaries and plot UMAP if one exists; it does not change fitted results. |
| `copy` | Return a modified copy if true; otherwise mutate input and return `None`. |

### Docstring outputs (1676-1704)

The function writes five columns to `adata.obs` and a fit dictionary to
`adata.uns['modality_em']`. With an `assay_col`, the labels describe agreement
or disagreement between assay origin and fitted component. Without it, the
labels are simply `Cell-like`, `Nucleus-like`, or `Neutral`.

The docstring also describes optional console summaries and UMAP plotting. A
missing UMAP does not prevent fitting.

### Validate thresholds and honor `copy` (1705-1711)

The only explicit ordering validation is
`0 <= gamma_lo < gamma_hi <= 1`. Strict separation prevents overlapping label
rules. The code does not similarly validate the anchor order, constraint order,
or `init` spelling.

### Build count evidence (1713-1723)

1. Sum raw sc-signature counts into `x` and sn-signature counts into `y`.
2. Warn if either contains non-integers, because that often means normalized or
   log-transformed input was supplied.
3. Round and cast both arrays to integers even after warning.
4. Compute `m=x+y` and `enough = (m >= min_counts)`.

Only `enough` observations fit the mixture and qualify for evidence-based final
labels. The function still computes posteriors for the other observations.

### Compress and fit globally (1725-1741)

- `uniq(mask)` compresses repeated `(x,m)` pairs into one row plus a frequency.
  This is computationally equivalent to fitting every repeated droplet but is
  much faster.
- Fewer than ten distinct informative pairs is an error.
- The selected initialization and constraints are passed to the global fit.
- If the one-component BIC is no worse than the two-component BIC, the function
  warns that the data may contain one modality. It intentionally continues and
  reports posteriors, so the user must inspect the fitted means and BIC.

### Optional cell-type refitting (1743-1779)

The global fit is always created first.

If `cell_type_col` is supplied, each cell type is handled as follows:

1. Require at least `min_cells_fit` informative observations.
2. Fit its unique `(x,m)` pairs, warm-starting from the global parameters.
3. Accept the local fit only if
   `mean(ab_cell) - mean(ab_nuc) >= min_separation`.
4. On too few observations, fitting failure, or insufficient separation, emit
   a warning and use the global fit.
5. Compute posteriors for all droplets in that cell type and tag them as
   `per_cell_type` or `global_fallback`.

Without `cell_type_col`, every posterior uses the global fit and
`modality_fit='global'`.

One subtlety: because local fits receive `init_from=fit_global`, the warm-start
branch takes precedence over a fresh anchored/split partition inside each cell
type. The chosen global initialization can still influence the local fit via
that starting point.

This is not a hierarchical or partially pooled statistical model. Each accepted
cell-type fit is optimized independently after the warm start. The global fit
provides initialization and an all-or-nothing fallback, but no penalty
continuously shrinks local parameters toward global parameters.

### Posterior thresholding (1781-1784)

```python
cell_p = gamma > gamma_hi
nuc_p  = gamma < gamma_lo
```

The comparisons are strict. Equality with a threshold remains undecided. Labels
start as `Neutral`, and subsequent rules overwrite only confident, informative
observations.

### Classification with assay metadata (1786-1819)

For recognized assay values, the mapping is:

| Known assay | Posterior condition | Final label |
|---|---|---|
| sc | `enough` and `p_cell > gamma_hi` | `Typical Cell (SC)` |
| sc | `enough` and `p_cell < gamma_lo` | `Nucleus-like Cell (SC)` |
| sn | `enough` and `p_cell < gamma_lo` | `Typical Nucleus (SN)` |
| sn | `enough` and `p_cell > gamma_hi` | `Cell-like Nucleus (SN)` |
| either in mixed data | insufficient evidence or middle posterior | `Neutral` |

If the recognized metadata contains only sc or only sn observations, lines
1799-1817 replace undecided labels belonging to that modality with its typical
label. Consequently, in recognized sc-only data, the anomalous label requires
a confident low posterior and enough counts; middle-posterior and low-count
droplets are folded into `Typical Cell (SC)`. The mirror rule applies to sn-only
data.

Values in `assay_col` that match neither `sc_label` nor `sn_label` remain
`Neutral`; line 1818 preserves that category if needed.

### Classification without assay metadata (1820-1823)

Confident informative droplets become `Cell-like` or `Nucleus-like`.
Everything else remains `Neutral`. There is no typical/anomalous interpretation
because no known assay identity was provided.

### Stored outputs and reporting (1825-1841)

| Location | Output | Interpretation |
|---|---|---|
| `obs['modality_counts']` | `m=x+y` | Amount of signature evidence. |
| `obs['sc_frac']` | `x/m`, or `NaN` when `m=0` | Raw cell-signature fraction. |
| `obs['p_cell']` | `gamma` | Model posterior for the high-mean component. |
| `obs['modality_fit']` | fit tag | Global, local, or fallback fit used. |
| `obs['modality_classification']` | categorical label | Thresholded result, optionally named using assay metadata. |
| `uns['modality_em']['global']` | fit dictionary | Global parameters, mixture weight, likelihood, and BIC. |
| `uns['modality_em'][cell_type]` | fit dictionary and tag | Present when per-cell-type fitting is requested. |

### How to interpret a stored fit

The raw dictionary is most readable after converting each `(a,b)` pair to mean,
concentration, and overdispersion correlation:

```python
fit = adata.uns["modality_em"]["global"]

a_cell, b_cell = fit["ab_cell"]
a_nuc, b_nuc = fit["ab_nuc"]

p_cell_component = a_cell / (a_cell + b_cell)
p_nuc_component = a_nuc / (a_nuc + b_nuc)

kappa_cell = a_cell + b_cell
kappa_nuc = a_nuc + b_nuc

rho_cell = 1 / (kappa_cell + 1)
rho_nuc = 1 / (kappa_nuc + 1)

separation = p_cell_component - p_nuc_component
delta_bic = fit["bic_1"] - fit["bic_2"]
```

Interpret the derived values as follows:

| Derived value | Statistical meaning |
|---|---|
| `p_cell_component` | Center of the fitted high-sc-fraction population. |
| `p_nuc_component` | Center of the fitted low-sc-fraction population. |
| `kappa_cell`, `kappa_nuc` | Within-component concentration; larger means less latent heterogeneity. |
| `rho_cell`, `rho_nuc` | Beta-binomial intra-class correlation/overdispersion; larger means more extra-binomial variation. |
| `separation` | Distance between component means on the sc-fraction scale. |
| `lam_cell` | Fitted prevalence of the high-mean component among informative observations. |
| `delta_bic` | Positive favors two components over one under the code's BIC calculation. |
| `loglik` | Maximized/last accepted two-component mixture log likelihood; useful mainly for comparing fits to the same data. |

Do not confuse `p_cell_component` with the per-droplet column `obs['p_cell']`.
The former is one location parameter for the entire component; the latter is a
different posterior probability for every droplet.

With `print_mode=True`, the function prints label counts, an assay-by-label
crosstab, and a classification UMAP when `adata.obsm['X_umap']` exists. It only
warns when the embedding is absent.

Finally, line 1841 follows the same mutation convention as the old function:

```python
result = Scanpyplus.CellorNucEM(..., copy=True)   # result is AnnData
Scanpyplus.CellorNucEM(adata, ..., copy=False)    # adata changes; return is None
```

## 10. Input contract for the supplied quiz files

Read-only schema inspection of the 6k files gives:

| Dataset | Shape | Signature overlap (sc/sn) | Count storage | Relevant metadata |
|---|---:|---:|---|---|
| mixed dataset A | 6,000 x 29,987 | 10/9 | integer `X` and integer `layers['counts']` | `suspension_type` contains `cell` and `nucleus`; `cell_type` is available |
| single-modality dataset B | 6,000 x 28,562 | 10/10 | integer `X` and integer `layers['counts']` | no pre-existing sc/sn modality column |

This has two API consequences:

1. `layer='counts'` is an explicit, readable way to satisfy the likelihood's
   raw-count assumption.
2. If dataset A's `suspension_type` is used as `assay_col`, its exact label
   values must be supplied as `sc_label='cell'` and `sn_label='nucleus'`.
   Dataset B can be run without `assay_col` for generic labels, or a known sc
   metadata column can be added as demonstrated in the official notebook when
   SC-specific label names are desired.

The `assay` column in dataset A records sequencing technologies (`10x 5' v1`
and `10x multiome`), not the exact default strings `sc` and `sn`; column meaning
matters more than the column name.

The mixed 6k file lacks one provisional sn gene, `AC098829.1`. This is allowed:
`_set_counts_EM` uses the nine sn genes that are present. The overlap count is a
useful item to record because a different gene index or annotation can change
how much of a supplied signature is measurable.

## 11. A disciplined reading/debugging checklist

Before interpreting any run, verify:

1. Both signature lists overlap `adata.var_names` and note how many genes were
   found; `_set_counts_EM` silently allows partial overlap.
2. The selected matrix/layer contains raw integer counts.
3. `assay_col`, `sc_label`, and `sn_label` agree exactly if metadata-aware names
   are requested.
4. How many droplets satisfy `modality_counts >= min_counts`.
5. The global component means
   `a_cell/(a_cell+b_cell)` and `a_nuc/(a_nuc+b_nuc)`.
6. Whether `bic_2 < bic_1`, while remembering this is model-selection evidence,
   not biological ground truth.
7. Whether a cell-type result used a local fit or `global_fallback`.
8. Whether a reported “typical” label in single-modality data was a confident
   high/low call or an undecided label reassigned by lines 1799-1817.
9. Whether `copy=True` was used before assigning the returned value.

## 12. Compact comparison

| Feature | `CellorNuc` | `CellorNucEM` |
|---|---|---|
| Core signal | Difference of Scanpy gene-set scores | Counts and fraction across the two signatures |
| Expected scale | Active expression matrix; not checked | Raw integer counts; warns and rounds otherwise |
| Uses assay labels to fit? | Yes | No |
| Uses cell types? | Required for thresholds | Optional local refits |
| Uncertainty | Between empirical score thresholds | Posterior between `gamma_lo` and `gamma_hi`, plus low evidence |
| Model diagnostics | None | Component parameters, likelihood, BIC, fit source |
| Quiz relevance | Background/legacy contrast | Primary implementation |
