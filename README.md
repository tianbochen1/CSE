# Copula-spectral estimation

Code and numerical results for **Copula-Spectral Estimation of Nonlinear Time-Series Models**, by Tianbo Chen and Ruidi Zhang.

The method estimates Gaussian stochastic-volatility dependence parameters by matching empirical integrated rank spectra to deterministic model spectra. The marginal distribution is unspecified. The parameters are latent persistence `phi` and stationary latent variance `v`.

## Contents

| Path | Purpose |
| --- | --- |
| `R/cse_core.R` | Hermite coefficients, empirical spectra, SV simulation, estimation, and comparison methods |
| `R/cse_experiments.R` | Estimation, marginal transformations, feature selection, and numerical refinement |
| `R/run_simulation.R` | Estimation experiments with 500 replications per setting |
| `R/cse_applications.R` | Price-data preparation, dense-grid estimation, and spectral checks |
| `R/run_applications.R` | EURUSD and S&P 500 fits, tie sensitivity, and refitted S&P 500 checks |
| `python/cse.py` | Model evaluation and refitted parametric bootstrap |
| `python/run_inference.py` | Interval coverage and model-check calibration study |
| `python/manuscript_outputs.py` | All four figures, interval summaries, and pointwise reference envelopes |
| `python/application_sensitivity.py` | Application grid refinements and EURUSD marginal-transform QML fits |
| `python/check_numerics.py` | Numerical precision checks at the application parameter bound |
| `results/reported/` | Numerical summaries and spectral curves |
| `results/inference_replicates.csv` | All 12,000 outcomes of the inference study |
| `results/inference_table.csv` | Generated coverage, mean interval length, rejection frequency, and boundary frequency |
| `results/SP500_reference_draws.csv` | Generated empirical lower-upper quantile spectra underlying the pointwise envelopes |
| `results/sensitivity/` | Generated application grid and marginal-transformation results |

## Estimation experiments

The R analysis uses base R and the recommended `parallel` package. Run from the repository root:

```sh
Rscript R/run_simulation.R
```

The main experiment uses 500 replications, sample sizes 256, 512, 1024, and 2048, and the region `[0.005, 0.98] x [0.01, 3]`. The parameter grid has persistence values 0.3, 0.6, and 0.9 and latent variances 0.25, 1, and 2. Comparisons within an experiment use the same simulated paths. Separate experiments examine increasing marginal transformations, quantile pairs, integration endpoints, criterion-grid refinement, and pairwise likelihood. Settings and worker counts are specified in the runner.

Results under `results/reported/` are numerical summaries, including rounded table entries. The baseline file covers all 36 settings. The feature and marginal files include sample sizes omitted from the printed tables. The estimation runner saves individual records and computes its summaries from them. Additional block-comparison and timing summaries are also retained in this directory.

## Data and applications

The data sources are FRED series [DEXUSEU](https://fred.stlouisfed.org/series/DEXUSEU), supplied by the Board of Governors of the Federal Reserve System, and [SP500](https://fred.stlouisfed.org/series/SP500), supplied by S&P Dow Jones Indices LLC. The exchange rate is in US dollars per euro; the equity series is a price index excluding dividends.

Returns are `100 * diff(log(price))` over successive available prices, dated 3 January 2017 through 31 December 2025, with the initial price dated 30 December 2016. The samples contain 2246 and 2262 returns. Raw provider prices are not redistributed here; the return series used in the study are available from the authors on request.

```sh
Rscript R/run_applications.R
python python/application_sensitivity.py --data-dir data
```

The R runner downloads the FRED series into `data`, or accepts complete local FRED-format CSV files through `LOCAL_CSV`. The sensitivity program takes `DEXUSEU.csv` and `SP500.csv` from its data directory. Historical provider coverage can change; local files must cover the stated dates.

Application fits use `[0.005, 0.995] x [0.01, 3]`, with 128 quantile points and 992 integration endpoints for both full samples. The sensitivity program reproduces six grid resolutions, evaluates their targets on the reference grid, and fits QML to the four centered EURUSD marginal transformations. Its criterion uses the same uniformly weighted spectrum as the R application estimator, evaluated through matrix products.

EURUSD is used for descriptive estimation and sensitivity to upper, average, and randomized tied ranks. S&P 500 model checks use 499 stationary samples, with ranks recomputed and parameters re-estimated in every sample. The global statistic is the sum of its real and imaginary components. The median comparison uses 20 indicator autocorrelations. The bootstrap probabilities use the plus-one correction.

## Inference study

Install the Python dependencies and run:

```sh
python -m pip install -r requirements.txt
python python/run_inference.py --replications 1000 --bootstrap 199 --workers 8
python python/manuscript_outputs.py intervals
```

There are 12 settings: sample sizes 512, 1024, and 2048 crossed with `(phi, v)` equal to `(0.6, 1)`, `(0.9, 1)`, `(0.98, 0.25)`, and `(0.98, 1)`. Each setting has 1000 independent outer replications and 199 refitted bootstrap samples per replication. Basic and percentile intervals have nominal coverage 95%; checks use level 5%. Basic intervals are intersected with the parameter range. Coverage and mean length are reported together.

The settings file contains the design and seed. Independent streams are indexed by setting and replication. Numerical convergence and boundary flags are retained with every outcome. Monte Carlo standard errors for proportions are `sqrt(p * (1 - p) / R)`.

## Figures and reference envelopes

```sh
python python/manuscript_outputs.py figures
python python/manuscript_outputs.py envelopes
python python/plot_applications.py
```

The first command constructs Figures 1--4 from model evaluations and the stored summaries, generating the reference envelopes before plotting. Figures 2 and 3 use the rounded published RMSE values. The second command computes only the S&P 500 reference envelopes and saves their empirical spectral draws; the third generates the envelopes and plots the application figure. `python python/manuscript_outputs.py all` runs the summaries, envelopes, and all figures in that order. Generated files are written locally; the reported numerical summaries are left unchanged.

The Figure 4 envelopes are pointwise 2.5% and 97.5% quantiles of empirical centered integrated spectra from 499 stationary samples generated at the full-sample S&P 500 estimate. They treat that estimate as a fixed reference; parameters are not refitted for these envelopes. Quantiles use linear interpolation. The envelope command writes the sample size, fitted parameter, seed, and construction to `results/reference_envelope_settings.json`, the 499 spectral curves to `results/SP500_reference_draws.csv`, and the plotted curves to `results/reference_curves/`. The model-check probabilities in the application table use parameter re-estimation and the full evaluation grid.

## Numerical accuracy

```sh
python python/check_numerics.py
```

This compares interpolation with direct coefficient evaluation at 17 variance values, all pairs of the 128 application quantiles, and five integration endpoints, using persistence 0.995. It also compares 1856 with 1984 Gaussian-quadrature points and Hermite orders 1791 with 2251. The resulting errors and exact settings are saved in `results/application_numerical_accuracy.csv` and `results/application_numerical_settings.json`.
