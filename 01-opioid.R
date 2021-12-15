### Opioid treatment data, town-level random intercept ###
# This script implements data cleaning and runs ELGM model fits
# for Opioid example
# This script may be sourced in a clean R session, after setting the paths
# found in the SETUP section.

## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'Matrix',
  'aghq',
  'TMB',
  'glmmTMB'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

TMB::precompile()

# Check that OpenMP is installed and working.
# The time comparisons in the paper require this. If not, the times will be slower
# than what we reported.
stopifnot(TMB::openmp() > 0)

## Global variables ----
# where to save the results
resultspath <- "~/data/"

# How many times to do each simulation
numtimes <- 10

# Set this flag if you want to process and save the data- default FALSE
# means this script just sources functions necessary for model fitting
dodata <- TRUE

## SET THESE PATHS ACCORDING TO WHERE THE DATA AND CODE ARE SAVED ##
# Full path of the data. Change this to wherever you want the data saved to
datapath <- "/home/alex/data/30122-0001-Data.tsv"

# Path where you want to save the data as RDS
savepath <- "/home/alex/data/opioid-admission-data-rds.rds"

# Compile and load the TMB template. Source code included in the aghq package
tmbpath <- tempdir()
stopifnot(file.copy(system.file('extsrc/01_opioid.cpp',package='aghq'),tmbpath))
compile(file.path(tmbpath,"01_opioid.cpp"))
dyn.load(dynlib(file.path(tmbpath,"01_opioid")))

## END SETUP ##

## Read in and preprocess data ----

if (dodata){
  ## Read in the data ----
  # This step only has to be done once, will save as rds for faster loading
  # during the time comparisons
  
  opioid <- readr::read_tsv(
    file = datapath,
    col_names = TRUE,
    col_types = stringr::str_c(rep("n",65),collapse = "")
  )
  # glimpse(opioid)
  
  ## Clean the data ----
  
  # Subset only variables we want, and create the response: ----
  # an indicator of whether the treatment was completed.
  # Also remove missing values
  # Also combine some groups and recode levels to set
  # sensible references
  opioid_clean <- opioid %>%
    dplyr::select(
      year = DISYR,
      id = CASEID,
      age = AGE,
      gender = GENDER,
      race = RACE,
      education = EDUC,
      employment = EMPLOY,
      state = STFIPS,
      town = CBSA,
      livingarrangement = LIVARAG,
      reason = REASON
    ) %>%
    mutate_at(c("age","gender","race","education","employment","state","town","livingarrangement"),as.character) %>%
    mutate(died = if_else(reason == 6,1,0),
           completed = if_else(reason == 1,1,0)) %>%
    filter(reason != -9,gender != "-9",race != "-9",state != "-9",livingarrangement != "-9",town != "-9") %>%
    mutate(
      gender = case_when(
        gender == "1" ~ "0M",
        gender == "2" ~ "F",
        TRUE ~ "-9"
      ),
      livingarrangement = case_when(
        livingarrangement == "1" ~ "H",
        livingarrangement == "2" ~ "D",
        livingarrangement == "3" ~ "0I",
        TRUE ~ "-9"
      ),
      race = case_when(
        race %in% c("1","2","23") ~ "AI",
        race %in% c("3","13") ~ "AS",
        race == "4" ~ "B",
        race == "5" ~ "0W",
        race %in% c("20","21") ~ "O",
        TRUE ~ "-9"
      )
    )
  
  # Save the data to rds for faster loading
  readr::write_rds(opioid_clean,savepath,compress = "none")
}
## Function to fit model ----

fit_model <- function(dat) {
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
  
  quadmod <- aghq::marginal_laplace_tmb(template,7,template$par,control = default_control_tmb(method='BFGS'))
  # Include time taken for summaries and posterior samples in the time comparison
  # aghq:::summary.marglaplace (with max_print = Inf) takes care of this
  summary(quadmod,max_print = Inf)
  
  quadmod
}

## BEGIN TIMING ##
# Perform timings: run ELGM numiter times for various subsamples of the opioid data
# Do not save the results to disk; creation of the model objects used for comparing
# outputs is done in 01-opioid-compare.R

cat("Reading in data, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

opioid_clean <- read_rds(savepath)

cat("Starting sims, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

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

tm_nfull <- numeric(numtimes)
for (tt in 1:numtimes) {
  tm <- Sys.time()
  model_obj_nfull <- fit_model(opioid_clean)
  tm_nfull[tt] <- difftime(Sys.time(),tm,units = "secs")
  cat("n = full, time = ",tm_nfull[tt],", iteration",tt,".\n")
}

cat("Done timing. Saving results.\n")
## END TIMING ##

timeframe <- tibble(
  itr = rep(1:numtimes,5),
  n = rep(c(10^(3:6),nrow(opioid_clean)),each=numtimes),
  time = c(tm_n1000,tm_n10000,tm_n100000,tm_n1000000,tm_nfull)
)

print(timeframe,n = Inf)

write_csv(timeframe,file.path(resultspath,"timeframe-20211022.csv"))
cat("Done.\n")
