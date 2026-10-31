# Copula-spectral estimation

Code and numerical summaries for **Copula-Spectral Estimation of Nonlinear Time-Series Models**, by Tianbo Chen and Ruidi Zhang.

The method estimates Gaussian stochastic-volatility dependence parameters by matching empirical rank-based integrated copula spectra to deterministic model spectra. The marginal distribution is unspecified. The two parameters are latent persistence `phi` and stationary latent variance `v`.

## Contents

| Path | Purpose |
| --- | --- |
| `R/cse_core.R` | Hermite coefficients, empirical spectra, SV simulation, minimum-distance and comparison estimators |
| `R/cse_experiments.R` | Estimation, marginal transformations, feature selection, quadrature, and timing experiments |
| `R/run_simulation.R` | Run the estimation experiments with 500 replications per setting |
| `R/cse_applications.R` | Price-data preparation, dense-grid estimation, and diagnostic functions |
| `R/run_applications.R` | EURUSD and S&P 500 estimation, rolling comparisons, tie sensitivity, and S&P 500 model checks |
| `python/cse.py` | Vectorized model evaluation and refitted parametric bootstrap |
| `python/run_inference.py` | Coverage and model-check calibration study |
| `python/plot_applications.py` | Recreate the application figure from the reported spectral curves |
| `results/reported/` | Numerical summaries underlying the estimation study and applications |
| `results/inference_replicates.csv` | Individual outcomes of the inference study |
| `results/inference_summary.csv` | Coverage, interval length, bias, standard errors, rejection rates, and numerical diagnostics |
| `results/inference_settings.json` | Inference-study design and random seed |

## R experiments

The R analysis uses base R and the recommended `parallel` package. Run from the repository root:

```sh
Rscript R/run_simulation.R
Rscript R/run_applications.R
```

The estimation experiment uses 500 replications, sample sizes 256, 512, 1024, and 2048, and the parameter rectangle `[0.005, 0.98] x [0.01, 3]`. The full comparison includes persistence values 0.3, 0.6, and 0.9 and latent variances 0.25, 1, and 2. Separate experiments evaluate increasing marginal transformations, quantile pairs, integration endpoints, grid refinement, and pairwise likelihood. Experiment codes are listed at the top of the runner; `E01` is population identification, `E02` the full parameter grid, `E03` marginal transformations, `E04` quantile selection, `E05` endpoint selection, `E10` quadrature sensitivity, and `E11` likelihood refinement and timing. Replications and worker counts are explicit settings at the top of each runner.

The reported estimation-study timings were obtained with R 4.2.3 on 64-bit Windows. Timing ratios depend on the implementation and hardware. The cleaned R modules also pass parsing and numerical smoke checks with R 4.6.0 under WebAssembly: the accelerated and explicit spectral criteria agree to below `5e-17` in the checked cases. This smoke check does not replace execution of the full application analysis.

## Data and applications

The sources are FRED series [DEXUSEU](https://fred.stlouisfed.org/series/DEXUSEU) (US dollars per euro; Board of Governors of the Federal Reserve System) and [SP500](https://fred.stlouisfed.org/series/SP500) (price index; S&P Dow Jones Indices LLC). Returns are `100 * diff(log(price))` over successive available prices, dated 3 January 2017 through 31 December 2025, with the initial price dated 30 December 2016. The reported samples contain 2246 and 2262 returns, respectively.

The application runner downloads the price series or accepts complete local FRED-format CSV files through `LOCAL_CSV`. It checks date coverage and missing observations. FRED access and historical coverage, especially for the S&P 500 series, may change. The exact return series used for the reported results are available from the authors on request. Raw provider price data are not redistributed here.

Application fits use the rectangle `[0.005, 0.995] x [0.01, 3]`, 16 times the default number of quantile points, and 32 times the default number of integration endpoints. The full-sample grid has 128 quantile points and 992 endpoints. Rolling evaluation uses four common blocks of 256 returns and training lengths 512 and 1024. All comparison methods and the independence benchmark are retained in the numerical summaries.

EURUSD is used for descriptive dependence estimation and sensitivity to upper, average, and randomized ranks. Its rounded observations contain ties; continuous-model bootstrap p-values and reference envelopes are not reported for this series. The S&P 500 model checks use 499 stationary SV samples, with ranks recomputed and parameters re-estimated for every sample. The diagnostics assess the fitted SV copula; the imaginary component is not a test against every possible reversible process.

## Python inference study

Use Python 3.11 or later and the versions in `requirements.txt`:

```sh
python -m pip install -r requirements.txt
python python/run_inference.py --replications 1000 --bootstrap 199 --workers 8
python python/plot_applications.py
```

The inference design contains 12 settings: sample sizes 512, 1024, and 2048 crossed with `(phi, v)` equal to `(0.6, 1)`, `(0.9, 1)`, `(0.98, 0.25)`, and `(0.98, 1)`. Each setting has 1000 independent outer replications and 199 bootstrap samples per replication. Every bootstrap sample is reranked and fully refitted. No linearized or one-step refit is used. The parameter rectangle is `[0.005, 0.995] x [0.01, 3]`. Basic and percentile intervals have nominal coverage 95%; the global, real-component, and imaginary-component checks use nominal level 5%.

The seed is fixed in the settings file. Independent streams are indexed by setting and replication, so the result does not depend on the order of worker completion. Completed replications can be resumed with the same settings. A different design requires a different output directory. The coefficient bank is generated automatically and is not versioned.

A 25 by 25 parameter grid initializes local optimization in the Python bootstrap implementation. In 36 sampled comparisons, it achieved the same minimum as the four-start estimator within `1.6e-17` in objective value. Analytic derivatives were checked against central differences; the maximum discrepancy was below `1.7e-12`. These checks concern numerical equivalence, not statistical coverage.

## Interpretation of stored results

Files under `results/reported/` contain published numerical summaries, including rounded table entries. They are not individual Monte Carlo records. The baseline file includes all 36 sample-size and parameter settings, including weak dependence. The feature, marginal, and timing files retain sample sizes omitted from the compact printed appendix.

The inference files contain all outer replications and both interval constructions, including boundary estimates and numerical convergence flags. Monte Carlo standard errors for a proportion `p` are `sqrt(p*(1-p)/R)`. Finite-sample interval performance varies substantially across settings, so interval coverage and model-check size should be assessed separately.

The application-curve files contain the empirical lower-upper quantile spectrum and the fitted model target. For the S&P 500 they also contain pointwise 2.5% and 97.5% reference quantiles. These envelopes are not simultaneous confidence bands.
