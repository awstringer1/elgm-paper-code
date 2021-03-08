# Code for Fast, Scalable Approximations to Posterior Distributions in Extended Latent Gaussian Models

This README describes the code files in this repository which replicate the examples from the manuscript *Fast, Scalable Approximations to Posterior Distributions in Extended Latent Gaussian Models* (arXiv).

All examples make use of the `aghq` package for approximate Bayesian inference using adaptive Gauss-Hermite quadrature, which is on `CRAN` ([github](https://github.com/awstringer1/aghq/), [arXiv](https://arxiv.org/abs/2101.04468)). See also this [paper](https://arxiv.org/abs/2102.06801) and associated [code](https://github.com/awstringer1/aghq-paper-code).

## Data

Data are sourced as follows:

- **Section 5**: data are sourced from the [ICPSR](https://www.icpsr.umich.edu/web/ICPSR/studies/30122), choose "Download" and then "Delimited". The user agreement prohibits directly sharing the data here.

- **Section 6.1**: data are available [here](https://github.com/aknandi/disaggregation_paper/tree/master/data).

- **Section 6.2**: dataset `Leuk` in package `INLA`.

- **Section 6.3**: data are printed in Table 1 of *BAYESIAN MASS ESTIMATES OF THE MILKY WAY: THE DARK AND LIGHT SIDES OF PARAMETER ASSUMPTIONS* by Eadie and Harris, *The Astrophysical Journal*, 829:108, and available as objects `gcdata` and `gcdatalist` in package `aghq`.

## Scripts

The examples are replicated in the following scripts. In all cases you should change the paths and flags at the top of the files as appropriate.

- **Section 5**:
  - Setup and data cleaning: *01-opioid-data-and-functions-inla.R*,
  - Run models and compare timings: *01-opioid-time-compare.R*.

- **Section 6.1**:
  - ELGM: *02-disaggregation.R*,
  - MCMC: *02-disaggregation-mcmc.R*.

- **Section 6.2**: *03-coxph.R*.

- **Section 6.3**:
  - Inference: *04-astro.R*,
  - `TMB` function template: `13_astro.cpp`.

Packages required: `tidyverse`, `disaggregation`, `Matrix`, `aghq`, `mapmisc`, `ipoptr`, `trustOptim`, `TMB`, `geostatsp`, `sp`, `raster`, `tmbstan`, `INLA`.
