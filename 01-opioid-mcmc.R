### Opioid treatment data, town-level random intercept ###
# This script fits the opioid model using tmbstan

## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'Matrix',
  'aghq',
  'TMB',
  'glmmTMB',
  'tmbstan'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

TMB::precompile()

## Global variables
# where to save the results
savepath <- "~/data/

# Number of total iterations. This need NOT be large, since the model objects created
# in this script are only ever used for time comparison, which is done on a per-iteration basis.
# A larger number of iterations are used for the output comparisons in 01-opioid-compare.R
numiter <- 1000

# Full path of the preprocessed data
datapath <- "/home/alex/data/opioid-admission-data-rds.rds"

# Compile and load the TMB template. Source code included in the aghq package
tmbpath <- tempdir()
stopifnot(file.copy(system.file('extsrc/01_opioid.cpp',package='aghq'),tmbpath))
compile(file.path(tmbpath,"01_opioid.cpp"))
dyn.load(dynlib(file.path(tmbpath,"01_opioid")))

## END SETUP ##

## Read in data ----

opioid_clean <- read_rds(datapath)

# Write a function to fit the model
fit_model <- function(dat,iter=numiter) {
  # iter: total number of iterations
  # will use ceiling(iter/10) as warmup
  
  template_info <- glmmTMB::glmmTMB(
    completed ~ gender + race + livingarrangement + (1|state) + (1|town),
    data = dat,
    family = binomial,
    doFit = FALSE,
    sparseX = c('cond' = TRUE)
  )

  template_data <- with(template_info$data.tmb,list(
    X = as(XS,'dgTMatrix'),
    Z = as(Z,'dgTMatrix'),
    y = yobs
  ))
  template_param <- list(
    beta = rep(0,ncol(template_data$X)),
    Us = rep(0,length(unique(dat$state))),
    Ut = rep(0,length(unique(dat$town))),
    theta = rep(0,2)
  )
  stopifnot(length(template_param$Us) + length(template_param$Ut) == ncol(template_data$Z))
  stopifnot(length(template_param$beta) == ncol(template_data$X))
  with(template_data,stopifnot(nrow(X) == nrow(Z)))
  with(template_data,stopifnot(nrow(X) == length(y)))
  
  template <- TMB::MakeADFun(
    data = template_data,
    parameters = template_param,
    random = c('beta','Us','Ut'),
    silent = TRUE,
    DLL = "01_opioid"
  )
  
  stanmod <- tmbstan::tmbstan(
    template,
    chains = 1,
    iter = iter,
    warmup = ceiling(iter/10)
  )
  # Convert to samples, as part of timing.
  as.data.frame(stanmod)
  # ...but actually return the whole object, for later processing
  stanmod
}

## BEGIN TIMING ##

cat("Starting timings, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

set.seed(8967687)
tm <- Sys.time()
model_obj_n1000 <- fit_model(sample_n(opioid_clean,1e03))
tm_n1000 <- difftime(Sys.time(),tm,units = "secs")
cat("n = 1,000, time = ",tm_n1000,"\n")

tm <- Sys.time()
model_obj_n10000 <- fit_model(sample_n(opioid_clean,1e04))
tm_n10000 <- difftime(Sys.time(),tm,units = "secs")
cat("n = 10,000, time = ",tm_n10000,"\n")

tm <- Sys.time()
model_obj_n100000 <- fit_model(sample_n(opioid_clean,1e05))
tm_n100000 <- difftime(Sys.time(),tm,units = "secs")
cat("n = 100,000, time = ",tm_n100000,"\n")

tm <- Sys.time()

model_obj_n1000000 <- fit_model(sample_n(opioid_clean,1e06))
tm_n1000000 <- difftime(Sys.time(),tm,units = "secs")
cat("n = 1,000,000, time = ",tm_n1000000,"\n")


## END TIMING ##
cat("Done timing. Saving results.\n")

timeframe <- tibble(
  n = 10^(3:6),
  time = c(tm_n1000,tm_n10000,tm_n100000,tm_n1000000)
)

print(timeframe,n = Inf)

write_csv(timeframe,file.path(savepath,"mcmc-timeframe-20211022.csv"))
cat("Done.\n")
