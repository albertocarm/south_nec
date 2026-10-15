# SOUTH-NEC

Bayesian promotion-time cure model of the GETNE SOUTH-NEC study.

> Carmona-Bayonas A, Jimenez-Fonseca P, Capdevila J, et al. Perioperative chemotherapy and long-term cure in resected grade 3 gastroenteropancreatic neuroendocrine carcinomas: the GETNE SOUTH-NEC study. *European Journal of Cancer*, 2026. https://www.sciencedirect.com/science/article/pii/S0959804926015558

## Model

- Population survival S(t) = exp{−θ F(t)}, with Weibull promotion time F; cure probability exp(−θ).
- Clonogenic component, log θ: TNM stage, colorectal primary, Ki-67 and age (restricted cubic splines, knots 30, 70, 90% and 49, 71, 83 years), perioperative chemotherapy.
- Kinetic component, log λ: Ki-67 (same spline), perioperative chemotherapy.
- Chemotherapy is a time-dependent switch at its start date: the burden of the foci still dormant is multiplied by exp(δ) and their activation scale by exp(η).
- Missing Ki-67 and age are bounded parameters of the model.
- Priors on standardized covariates: Normal(0, 2.5) clonogenic, Normal(0, 0.25) kinetic, Normal(0, 2.5) intercepts; Gamma(2, 1) Weibull shape.
- Stan: 4 chains, 4,000 iterations (2,000 warmup), adapt_delta 0.95, max_treedepth 12, seed 42.

## Run

1. Download the repository (Code > Download ZIP) and unzip it.
2. Open R in that folder, or `setwd()` to it.
3. Run:

```r
source("R/cure_model.R")   # Table 4, Supplementary Tables A5 and A6, Figure 2, Supplementary Figure A3
source("R/sensitivity.R")  # Supplementary Tables A1 and A10, truncated follow-up
```

Data are read from `data/southnec.csv`; results are written to `results/`.

Requires R ≥ 4.1 with rstan, survival, loo, ggplot2 and patchwork, and a C++ toolchain for Stan (Rtools on Windows).

## Data

`data/southnec.csv` is a simulated toy dataset, provided only to run the code. It contains no patient data and does not reproduce the results of the article. The study data are available from GETNE upon request to the corresponding authors.

| Variable | Description |
|---|---|
| `time` | Months from surgery to the first event or last follow-up |
| `efs` | 1 = recurrence or death, 0 = censored |
| `rfs` | 1 = recurrence, 0 = censored |
| `chemo` | 1 = perioperative chemotherapy |
| `chemo_type` | none, adjuvant, neoadjuvant |
| `chemo_start` | Months from surgery to the start of chemotherapy (≤ 0 if before surgery); empty if none or after the first event |
| `stage` | TNM stage: I, II, III |
| `colorectal` | 1 = colorectal primary |
| `ki67` | Ki-67 index (%); empty if not recorded |
| `ki67_min` | Lower bound of a missing Ki-67 (%) |
| `age` | Age (years); empty if not recorded |

## Files

| File | Content |
|---|---|
| `stan/cure_model.stan` | Model |
| `R/setup.R` | Packages, data, Stan data, posterior summaries |
| `R/cure_model.R` | EFS and RFS models, diagnostics, posterior predictive check, Ki-67 × chemotherapy |
| `R/sensitivity.R` | Priors, linear Ki-67 or age, truncated follow-up |

## License

MIT
