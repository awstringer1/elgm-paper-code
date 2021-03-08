### Time comparisons for Opioid dataset with INLA and ELGM ###

# Set the libPaths to my main user
# Had to use an alternate user to get PARDISO to work
.libPaths("/home/alex/R/x86_64-pc-linux-gnu-library/4.0")

library(tidyverse)
library(Matrix)
library(aghq)
library(INLA) # Testing version installed on 2021/02/01
inla.setOption(pardiso.license="/home/forpardiso/sys/licenses/pardiso.lic")
stopifnot(inla.pardiso.check() == "SUCCESS: PARDISO IS INSTALLED AND WORKING")

## Global variables ----

# Which comparisons to do?
doelgm <- TRUE
doinlamulti <- TRUE
doinlasingle <- TRUE
docompare <- TRUE

# How many threads for INLA multi?
numthreads <- 16 # Used for both num.threads and blas.num.threads

# Set random seed
seed <- 8967687

# Set this flag if you want to process and save the data- default FALSE
# means that data are loaded from pre-saved RDS object, faster
dodata <- FALSE

# Full path of the data. Change this to wherever you want the data saved to
datapath <- "/storage/opioid/opioid-admission-data-2006-2011.tsv"

# Path where you want to save/load the data as RDS
savepath <- "/storage/opioid/opioid-admission-data-rds.rds"

# Output path for timings and data
outputpath <- "/storage/forpardiso/"

# Datestamp
datestamp <- "20210212"

## Source function script ----
source("/home/forpardiso/01-opioid-data-and-functions-inla.R")


## Run comparisons! ----

cat("Starting sims, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

cat("Reading in data, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")
opioid_clean <- read_rds(savepath)

if (doelgm) {
  cat("Starting ELGM timings, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n1000 <- fit_model(sample_n(opioid_clean,1e03))
  tm_n1000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000, time = ",tm_n1000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n10000 <- fit_model(sample_n(opioid_clean,1e04))
  tm_n10000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 10,000, time = ",tm_n10000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n100000 <- fit_model(sample_n(opioid_clean,1e05))
  tm_n100000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 100,000, time = ",tm_n100000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n1000000 <- fit_model(sample_n(opioid_clean,1e06))
  tm_n1000000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000,000, time = ",tm_n1000000,"\n")

  tm <- Sys.time()
  set.seed(seed)
  model_obj_nfull <- fit_model(opioid_clean)
  tm_nfull <- difftime(Sys.time(),tm,units = "secs")
  cat("n = full, time = ",tm_nfull,"\n")
  
  timeframe <- tibble(
    n = c(10^(3:6),nrow(opioid_clean)),
    # n = 10^(3:5),
    time = c(tm_n1000,tm_n10000,tm_n100000,tm_n1000000,tm_nfull)
  )
  
  print(timeframe,n = Inf)
  
  write_csv(timeframe,paste0(outputpath,"ELGMtimings-",datestamp,".csv"))
  # Save the model objects too
  options(scipen = 10)
  save(list = paste0("model_obj_n",c(10^(3:6),"full")),file = paste0(outputpath,"ELGMmodels-",datestamp,".RData"))
  # save(list = paste0("model_obj_n",10^(3:5)),file = paste0(outputpath,"ELGMmodels-",datestamp,".RData"))
  
}

if (doinlamulti) {
  cat("Starting INLA multi timings, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n1000 <- fit_inla_multi(sample_n(opioid_clean,1e03),numthreads)
  tm_n1000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000, time = ",tm_n1000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n10000 <- fit_inla_multi(sample_n(opioid_clean,1e04),numthreads)
  tm_n10000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 10,000, time = ",tm_n10000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n100000 <- fit_inla_multi(sample_n(opioid_clean,1e05),numthreads)
  tm_n100000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 100,000, time = ",tm_n100000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n1000000 <- fit_inla_multi(sample_n(opioid_clean,1e06),numthreads)
  tm_n1000000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000,000, time = ",tm_n1000000,"\n")

  tm <- Sys.time()
  set.seed(seed)
  model_obj_nfull <- fit_inla_multi(opioid_clean,numthreads)
  tm_nfull <- difftime(Sys.time(),tm,units = "secs")
  cat("n = full, time = ",tm_nfull,"\n")
  
  timeframe <- tibble(
    n = c(10^(3:6),nrow(opioid_clean)),
    # n = 10^(3:5),
    time = c(tm_n1000,tm_n10000,tm_n100000,tm_n1000000,tm_nfull)
  )
  
  print(timeframe,n = Inf)
  
  write_csv(timeframe,paste0(outputpath,"INLAmultitimings-",datestamp,"-t",numthreads,".csv"))
  # Save the model objects too
  options(scipen = 10)
  save(list = paste0("model_obj_n",c(10^(3:6),"full")),file = paste0(outputpath,"INLAmultimodels-",datestamp,"-t",numthreads,".RData"))
  # save(list = paste0("model_obj_n",10^(3:5)),file = paste0(outputpath,"INLAmultimodels-",datestamp,"-t",numthreads,".RData"))
  
}

if (doinlasingle) {
  set.seed(seed)
  cat("Starting INLA single timings, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n1000 <- fit_inla_single(sample_n(opioid_clean,1e03))
  tm_n1000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000, time = ",tm_n1000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n10000 <- fit_inla_single(sample_n(opioid_clean,1e04))
  tm_n10000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 10,000, time = ",tm_n10000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n100000 <- fit_inla_single(sample_n(opioid_clean,1e05))
  tm_n100000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 100,000, time = ",tm_n100000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_n1000000 <- fit_inla_single(sample_n(opioid_clean,1e06))
  tm_n1000000 <- difftime(Sys.time(),tm,units = "secs")
  cat("n = 1,000,000, time = ",tm_n1000000,"\n")
  
  tm <- Sys.time()
  set.seed(seed)
  model_obj_nfull <- fit_inla_single(opioid_clean)
  tm_nfull <- difftime(Sys.time(),tm,units = "secs")
  cat("n = full, time = ",tm_nfull,"\n")
  
  timeframe <- tibble(
    n = c(10^(3:6),nrow(opioid_clean)),
    time = c(tm_n1000,tm_n10000,tm_n100000,tm_n1000000,tm_nfull)
  )
  
  print(timeframe,n = Inf)
  
  write_csv(timeframe,paste0(outputpath,"INLAsingletimings-",datestamp,".csv"))
  # Save the model objects too
  options(scipen = 10)
  save(list = paste0("model_obj_n",c(10^(3:6),"full")),file = paste0(outputpath,"INLAsinglemodels-",datestamp,"-t",numthreads,".RData"))
  # save(list = paste0("model_obj_n",10^(3:5)),file = paste0(outputpath,"INLAsinglemodels-",datestamp,"-t",numthreads,".RData"))
  
}

if (docompare) {
  cat("Comparing numerical results for ELGM and INLA.\n")
  # Load INLA and ELGM and compare the posterior means
  compareenvELGM <- new.env()
  compareenvINLA <- new.env()
  load(paste0(outputpath,"INLAmultimodels-",datestamp,"-t",numthreads,".RData"),env = compareenvINLA)
  load(paste0(outputpath,"ELGMmodels-",datestamp,".RData"),env = compareenvELGM)
  
  # ELGM
  sampsize <- 1e03
  set.seed(475687)
  postmeanELGM1000 <- apply(sample_marginal(compareenvELGM$model_obj_n1000,sampsize)$samps,1,mean)
  postmeanELGM10000 <- apply(sample_marginal(compareenvELGM$model_obj_n10000,sampsize)$samps,1,mean)
  postmeanELGM100000 <- apply(sample_marginal(compareenvELGM$model_obj_n100000,sampsize)$samps,1,mean)
  postmeanELGM1000000 <- apply(sample_marginal(compareenvELGM$model_obj_n1000000,sampsize)$samps,1,mean)
  postmeanELGMfull <- apply(sample_marginal(compareenvELGM$model_obj_nfull,sampsize)$samps,1,mean)
  
  # INLA
  getpostmeanINLA <- function(mod) {
    c(
      mod$summary.random$state$mean,
      mod$summary.random$town$mean,
      mod$summary.fixed$mean
    )
  }
  
  postmeanINLA1000 <- getpostmeanINLA(compareenvINLA$model_obj_n1000)
  postmeanINLA10000 <- getpostmeanINLA(compareenvINLA$model_obj_n10000)
  postmeanINLA100000 <- getpostmeanINLA(compareenvINLA$model_obj_n100000)
  postmeanINLA1000000 <- getpostmeanINLA(compareenvINLA$model_obj_n1000000)
  postmeanINLAfull <- getpostmeanINLA(compareenvINLA$model_obj_nfull)

  # Scaled L2 norm
  l2norm <- function(x,y) {
    stopifnot(length(x) == length(y))
    m <- length(x)
    sqrt(sum( (x-y)^2 )/m)
  }
  normtable <- tibble(
    n = c(10^(3:6),7283575),
    # n = 10^(3:5),
    L2 = c(
      l2norm(postmeanELGM1000,postmeanINLA1000),
      l2norm(postmeanELGM10000,postmeanINLA10000),
      l2norm(postmeanELGM100000,postmeanINLA100000),
      l2norm(postmeanELGM1000000,postmeanINLA1000000),
      l2norm(postmeanELGMfull,postmeanINLAfull)
    )
  )
  
  readr::write_csv(normtable,file = paste0(outputpath,"postmeancompare-",datestamp,"-t",numthreads,".csv"))
  cat("Done comparing.\n")
}
