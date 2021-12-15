### Opioid treatment data, town-level random intercept ###
# This script fits the opioid model using R-INLA
# You must have R-INLA installed, and PARDISO working

## BEGIN SETUP ##

## Libraries ----
# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'INLA'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    if (pkg == "INLA") {
      install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
      require(pkg,character.only = TRUE,quietly = TRUE)
    } else {
      install.packages(pkg)
      require(pkg,character.only = TRUE,quietly = TRUE)
    }
  }
}
# Note: the INLA developers strongly suggest using PARDISO for speed. This requires
# a license. Instructions available by typing INLA::inla.pardiso() (opens a web browser)
# PARDISO is used in the timings in the paper, so it is REQUIRED to run this script,
# else the timings will probably be much longer than what we reported.
# Set your path to the PARDISO license file here:
inla.setOption(pardiso.license="/home/alex/sys/licenses/pardiso.lic")
# Check if PARDISO is working properly.
stopifnot(inla.pardiso.check() == "SUCCESS: PARDISO IS INSTALLED AND WORKING")



## Global variables
# where to save the data
savepath <- "~/data/

# Number of threads to run INLA with. Use the same number as for running ELGM
threads <- 16
# Number of times to run each simulation
numtimes <- 10

# Full path of the preprocessed data
# Change this to whereever you saved the preprocessed data in 01-opioid.R
datapath <- "/home/alex/data/opioid-admission-data-rds.rds"

## END SETUP ##

## Read in data ----
cat("Reading in data.\n")

opioid_clean <- read_rds(datapath)

## Function to fit model ----
# Write a function to fit the model
fit_model <- function(dat,tr=threads) {
  # tr: number of threads to use for num.threads and num.blas.threads. If INLA
  # crashes, lower this number.
  inlamod <- tryCatch(inla(
      completed ~ gender + race + livingarrangement + 
        f(state,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))) +
        f(town,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))),
      data = dat,
      family = 'binomial',
      control.compute = list(
        openmp.strategy = 'pardiso',
        smtp = 'pardiso'
      ),
      control.inla = list(
        strategy = 'gaussian',
        int.strategy = 'ccd'
      ),
      num.threads = tr,
      blas.num.threads = tr
    ),
    error = function(e) cat(paste0("Error when fitting INLA multi threaded, n = ",nrow(dat),".\nThe error was: ",e,"\n"))
  )
  inlamod
}

## BEGIN TIMING ##
# Perform timings: run ELGM numiter times for various subsamples of the opioid data
# Do not save the results to disk; creation of the model objects used for comparing
# outputs is done in 01-opioid-compare.R

cat("Starting timings, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

set.seed(456132)
tm_n1000 <- numeric(numtimes)
for (tt in 1:numtimes) {
  tm <- Sys.time()
  model_obj_n1000 <- fit_model(sample_n(opioid_clean,1e03))
  tm_n1000[tt] <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000, time = ",tm_n1000[tt],", iteration",tt,".\n")
}

tm_n10000 <- numeric(numtimes)
for (tt in 1:numtimes) {
  tm <- Sys.time()
  model_obj_n10000 <- fit_model(sample_n(opioid_clean,1e04))
  tm_n10000[tt] <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 10,000, time = ",tm_n10000[tt],", iteration",tt,".\n")
}

tm_n100000 <- numeric(numtimes)
for (tt in 1:numtimes) {
  tm <- Sys.time()
  model_obj_n100000 <- fit_model(sample_n(opioid_clean,1e05))
  tm_n100000[tt] <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 100,000, time = ",tm_n100000[tt],", iteration",tt,"\n")
}

tm_n1000000 <- numeric(numtimes)
for (tt in 1:numtimes) {
  tm <- Sys.time()
  model_obj_n1000000 <- fit_model(sample_n(opioid_clean,1e06))
  tm_n1000000[tt] <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000,000, time = ",tm_n1000000[tt],",iteration",tt,".\n")
}

# Only run INLA on full data once, crashes with too many threads, and takes a long time
cat("Fitting INLA to full data now.\n")
tm <- Sys.time()
model_obj_nfull <- fit_model(opioid_clean,tr=1) # Run with 1 thread, crashes with 16
tm_nfull <- difftime(Sys.time(),tm,units = "secs")
cat("n = full, time = ",tm_nfull,".\n")

## DONE TIMING ##
cat("Done timing. Saving results.\n")

timeframe <- tibble(
  itr = c(rep(1:numtimes,4),1),
  n = c(rep(c(10^(3:6)),each=numtimes),nrow(opioid_clean)),
  time = c(tm_n1000,tm_n10000,tm_n100000,tm_n1000000,tm_nfull)
)

print(timeframe,n = Inf)

write_csv(timeframe,file.path(savepath,"inla-timeframe-20211022.csv"))
cat("Done.\n")
