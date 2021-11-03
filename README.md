# Code for Fast, Scalable Approximations to Posterior Distributions in Extended Latent Gaussian Models

This README describes the code files in this repository which replicate the examples from the manuscript *Fast, Scalable Approximations to Posterior Distributions in Extended Latent Gaussian Models* (arXiv).

All examples make use of the `aghq` package for approximate Bayesian inference using adaptive Gauss-Hermite quadrature, which is on `CRAN` ([github](https://github.com/awstringer1/aghq/), [arXiv](https://arxiv.org/abs/2101.04468)).

## Data

Data are sourced as follows:

- **Section 5**: data are sourced from the [ICPSR](https://www.icpsr.umich.edu/web/ICPSR/studies/30122), choose "Download" and then "Delimited". The user agreement prohibits directly sharing the data here.

- **Section 6.1**: data are available [here](https://github.com/aknandi/disaggregation_paper/tree/master/data), but they are downloaded and processed within the script so no need to do anything.

- **Section 6.2**: dataset `Leuk` in package `INLA`. Loaded in the script, no need to do anything.

- **Section 6.3**: data are printed in Table 1 of *BAYESIAN MASS ESTIMATES OF THE MILKY WAY: THE DARK AND LIGHT SIDES OF PARAMETER ASSUMPTIONS* by Eadie and Harris, *The Astrophysical Journal*, 829:108, and available as objects `gcdata` and `gcdatalist` in package `aghq`. Loaded in the script, no need to do anything.

## Scripts

The examples are replicated in the following scripts. Every effort has been made to make these as portable as possible. Specific points to pay attention to:

1. In all scripts for Section 5, you have to set the path to the data you downloaded, as well as paths to save the results. These are located between the `## BEGIN SETUP ##` and `## END SETUP ##` comments at the top of each script.
2. In the script for Section 6.1, data are downloaded from the internet within the script. If those data change location or are otherwise unobtainable at the time of running, then the script will not run. This is checked for within the script and informative errors are printed.
3. The script for Section 6.3 depends on the `IPOPT` software, which can be onerous to install.
4. The time comparisons assume you are working on an OpenMP-enabled machine. If you are not, the speeds for the ELGM will be slower (I assume).

Packages are loaded, and installed if not found, at the top of each script. Two non-standard installations:

1. Package `INLA` is not on `CRAN`. It is installed within the scripts. You need to use the `PARDISO` sparse matrix library if attempting to replicate the time comparisons in Section 5. Instructions for enabling this are found in the script where it is needed.
2. Package `ipoptr` is not on `CRAN`, and requires installation of the `IPOPT` optimization software. See https://coin-or.github.io/Ipopt/INSTALL.html.

- **Section 5**: run the following scripts in the order provided:
  - Data cleaning and ELGM timings: *01-opioid.R*,
  - INLA timings: *01-opioid-inla.R*,
  - MCMC timings: *01-opioid-mcmc.R*,
  - Output comparisons and creation of Tables 1 and 2: *01-opioid-combine.R*.

- **Section 6.1**:
  - ELGM: *02-disaggregation.R*,
  - MCMC: *02-disaggregation-mcmc.R*.

- **Section 6.2**: *03-coxph.R*.

- **Section 6.3**: *04-astro.R*
