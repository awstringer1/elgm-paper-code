### Opioid treatment data, town-level random intercept ###
# This script attempts to fit the opioid GLMM example using VI. Users can change the
# sample size, priors, and optimization parameters, to observe the challenges in
# doing this.


## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'Matrix',
  'aghq',
  'rstan',
  'rstanarm'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

## Global variables

# Where the opioid data are saved
datapath <- "/home/alex/data/opioid-admission-data-rds.rds"

opioid_clean <- read_rds(datapath)

# Set control parameters for the VI algorithm
controlparams <- list(
  # meanfield = fully factored Gaussian; fullrank = completely unspecified covariance
  algorithm = 'meanfield',
  # Standard deviation for the Gaussian prior on the intercept and other betas
  prior_sd_beta = 1,
  # Parameters for the prior on the covariance matrix, see https://mc-stan.org/rstanarm/articles/priors.html#how-to-specify-flat-priors-and-why-you-typically-shouldnt
  prior_covariance = c(1,1,1,-log(.5)/.5),
  iter = 1e04
)

# Fit to a subsample of size n
# Set n = 0 to use the full dataset
n <- 1e04
## END SETUP ##

# Subsample the data

if (n == 0) {
  dat <- opioid_clean
} else {
  dat <- sample_n(opioid_clean,n)
}

# Do the VI
stanvb <- rstanarm::stan_glmer(
  completed ~ gender + race + livingarrangement + (1|state) + (1|town),
  data = dat,
  family = binomial,
  prior = normal(0,controlparams$prior_sd_beta),
  prior_intercept = normal(0,controlparams$prior_sd_beta),
  prior_covariance = decov(regularization = controlparams$prior_covariance[1],concentration = controlparams$prior_covariance[2],shape = controlparams$prior_covariance[3],scale = controlparams$prior_covariance[4]),
  algorithm = controlparams$algorithm,
  sparse = TRUE,
  init = "0",
  iter = controlparams$iter
)
# Algorithm crashes, or runs to completion and produces Pareto k-hat >> 1
# Under the control params set above, algorithm converged in 6700 iterations with k-hat = 1.76
# An expert VI user could probably tune parameters to get better results here.