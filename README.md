# Rej-Detect (mtDNA Mixture Detection)

Rej-Detect is a prototype Bayesian framework for estimating the fraction of donor mitochondrial DNA (sample **D**) present in a mixed sample (**RDM**) that is dominated by recipient DNA (sample **R**). The intended use case is **very low donor fractions** (e.g., <2%).

The current implementation focuses on:

- Selecting **donor-informative mtDNA positions** from variant caller outputs (currently: `mutserve`-style TSV tables).
- Modeling read counts at informative sites with a **binomial likelihood**.
- Incorporating a **sequencing error rate** `ε` into the expected donor-allele observation probability.
- Estimating the donor mixture proportion `p` by **grid-based Bayesian inference** over a restricted range (default 0%–5%).

## Mathematical model (current)

For each informative position `i`:

- `k_i` = number of reads supporting the D-specific allele
- `n_i` = total reads at that position
- `p` = donor mixture proportion (unknown)
- `ε` = sequencing error rate (assumed known/fixed)

Expected frequency of observing the D-specific allele in the mixture:

`f(p) = p*(1-ε) + (1-p)*ε`

Likelihood at each position:

`k_i ~ Binomial(n_i, f(p))`

Assuming independence across informative sites, the total likelihood is the product across positions (implemented as a sum of log-likelihoods).

The posterior over `p` is computed on a grid with a **uniform prior** over `p ∈ [0, p_max]` (default `p_max = 0.05`).

## What is implemented today

### Refactored Python modules
Reusable code lives in:

- `src/mixture_bayes/`
  - `inference.py`
    - `log_likelihood(p, k_reads, n_reads, epsilon)`
    - `bayesian_estimation(k_reads, n_reads, epsilon, p_max=0.05, p_grid_size=4000, ci=0.95)`
  - `simulation.py`
    - `generate_informative_positions(genome_length, n_positions)`
    - `simulate_reads(p_true, n_positions, mean_depth, epsilon, ...)`
    - `simulate_and_estimate(...)`
  - `preprocessing.py`
    - `preprocess_variant_tables(recipient_file, donor_file, mixed_file, ...)`
      - Selects donor-informative positions using variant-level thresholds.
      - Converts RDM VAFs into `(k_reads, n_reads)`.
      - Uses per-row `Coverage` as `n_reads` if available; otherwise falls back to `mean_depth`.
  - `plotting.py`
    - `plot_posterior_distributions(results_list, ...)`
    - `plot_true_vs_estimated(results_list, ...)`

### Notebook entry point
- `src/Bayesian_Mixture_Estimation_Pipeline_Refactored.ipynb`
  - Simulation validation (0.1%, 0.5%, 1%, 2%)
  - Depth sensitivity analysis
  - Example preprocessing of `mutserve` tables (R, D, RDM) and running Bayesian inference

## Inputs / expected data

### Variant tables
The preprocessing function expects TSV files that include at least:

- `Pos`
- `Ref`
- `Variant`
- `VariantLevel`

If present, `Coverage` will be used as the per-site read depth for RDM.

## How to run

### Option A: Run the refactored notebook
1. Open `src/Bayesian_Mixture_Estimation_Pipeline_Refactored.ipynb`.
2. Update the file paths for `recipient_file`, `donor_file`, and `mixed_file`.
3. Run all cells.

### Option B: Use modules from Python
From a working directory that can import from `src/`:

```python
from mixture_bayes.preprocessing import preprocess_variant_tables
from mixture_bayes.inference import bayesian_estimation

informative_df, common_variants_df = preprocess_variant_tables(
    recipient_file="/path/to/R.tsv",
    donor_file="/path/to/D.tsv",
    mixed_file="/path/to/RDM.tsv",
    mean_depth=800,
    donor_homoplasmy_threshold=0.99,
    recipient_low_threshold=0.01,
)

res = bayesian_estimation(
    k_reads=informative_df["k_reads"].values,
    n_reads=informative_df["n_reads"].values,
    epsilon=0.02,
)

print(res["p_est_mean"], res["lower_bound"], res["upper_bound"])
```

## Current limitations

- Informative-site selection is **threshold-based** and does not yet implement a full “major allele comparison across all positions” rule.
- Heteroplasmy filtering is only partially captured by the donor homoplasmy threshold and recipient low threshold.
- Depth/QC filtering is minimal (no explicit filtering on low depth in R/D/RDM beyond what your thresholds imply).
- Sequencing error `ε` is **global constant**; per-site error `ε_i` is not implemented.
- Simulation currently simulates read counts at informative positions (not a full read-level mixture simulation).
- No automated benchmarking against other methods is included yet.

## Next version (planned features)

### More rigorous informative-site selection
- Explicitly compare **major alleles** between R and D at each position.
- Exclude positions with evidence of **heteroplasmy** in either sample.
- Enforce minimum depth thresholds for R, D, and RDM and/or filter by base quality / mapping quality.
- Optional site blacklists (low complexity regions, known problematic mtDNA sites, NUMTs-prone loci).

### Improved error model
- Support **per-site** error rates `ε_i` (or error rates by base / context).
- Incorporate strand bias and quality-aware error adjustments (where available).

### Evaluation and detection limit reporting
- Automated simulation sweeps over:
  - depth
  - error rate
  - number of informative positions
- Compute standardized metrics:
  - calibration / bias
  - sensitivity / specificity
  - minimum detectable mixture ratio under defined confidence criteria

### Pipelineization & reproducibility
- Batch-mode runner (CLI) to process many triplets (R, D, RDM).
- Standard output artifacts:
  - per-sample posterior plots
  - tables of estimates + credible intervals
  - QC summaries of informative sites used
- Add tests for core inference and preprocessing.

---

## Status
This repository currently contains a **working prototype** of the Bayesian estimator and a refactored notebook-driven workflow. It is ready for iterative refinement toward the full Rej-Detect specification.
