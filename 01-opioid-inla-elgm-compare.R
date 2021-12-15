### Combine results from all analyses of the Opioid data ###


## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'Matrix',
  'aghq',
  'TMB',
  'glmmTMB',
  'INLA',
  'tmbstan',
  'rstan'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    if (pkg == "INLA") {
      install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
    } else {
      install.packages(pkg)
      require(pkg,character.only = TRUE,quietly = TRUE)
    }
  }
}

TMB::precompile()
# Note: the INLA developers strongly suggest using PARDISO for speed. This requires
# a license. Instructions available by typing INLA::inla.pardiso() (opens a web browser)
# PARDISO is used in the timings in the paper.
# Comment the next two lines if you don't want to use PARDISO when running this code:
inla.setOption(pardiso.license="/home/alex/sys/licenses/pardiso.lic")
stopifnot(inla.pardiso.check() == "SUCCESS: PARDISO IS INSTALLED AND WORKING")


# The number of iterations used in the MCMC time comparisons
# This is only for calculating the "number of equivalent iterations" below.
numitr <- 1000 

## Set paths ----

globalpath <- "/home/alex/data"
datapath <- file.path(globalpath,"opioid-admission-data-rds.rds")

# Compile and load the TMB template. Source code included in the aghq package
tmbpath <- tempdir()
stopifnot(file.copy(system.file('extsrc/01_opioid.cpp',package='aghq'),tmbpath))
compile(file.path(tmbpath,"01_opioid.cpp"))
dyn.load(dynlib(file.path(tmbpath,"01_opioid")))

## END SETUP ##

## Fit all three methods to the same datasets ----
fit_all_models <- function(dat,domcmc=TRUE,tr=16,inlaintstrategy = "ccd",iter=10000) {
  # domcmc: set to FALSE to skip mcmc
  # tr: number of threads to use for INLA, set lower if it crashes
  ## ELGM ##
  cat("Fitting ELGM.\n")
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
  ## END ELGM ##
  
  ## INLA ##
  cat("Fitting INLA.\n")
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
      int.strategy = inlaintstrategy
    ),
    num.threads = tr,
    blas.num.threads = tr
  ),
  error = function(e) cat(paste0("Error when fitting INLA multi threaded, n = ",nrow(dat),".\nThe error was: ",e,"\n"))
  )
  ## END INLA ##
  
  ## MCMC ##
  if (domcmc) {
    cat("Fitting MCMC.\n")
    stanmod <- tmbstan::tmbstan(
      template,
      chains = 1,
      iter = iter,
      warmup = ceiling(iter/10)
    )
  } else {
    stanmod <- list()
  }
  ## END MCMC ##
  cat("Done.")
  list(
    aghq = quadmod,
    inla = inlamod,
    mcmc = stanmod
  )
}

get_samples <- function(mod,M=1e04) {
  # Compute summaries for each model
  # Return marginal posterior samples in a common format, for further downstream analysis
  # M: number of posterior samples to use for INLA and aghq
  # Note: "theta" means different things for each model because of internal parametrizations,
  # and the way different softwares pre-process the output. This function outputs
  # "theta" = "sigma", the standard deviation of the random effects (one for state,
  # and one for town)
  # Output: list(betasamps,thetasamps); matrices with samples in rows, parameters in columns
  
  # Switch based on which model
  if (inherits(mod,'aghq')) {
    ## W ##
    # get samples from JOINT posterior
    # Sorry for misleading function name-- see ?aghq::sample_marginal
    samps <- aghq::sample_marginal(mod,M = M)
    # Beta
    betasamps <- t(samps$samps[rownames(samps$samps) == 'beta', ])
    colnames(betasamps) <- paste0('beta',1:ncol(betasamps))
    # theta
    thetasamps <- samps$thetasamples
    # Transform from LOG-PRECISION to STANDARD DEVIATION
    thetasamps <- Reduce(cbind,lapply(thetasamps,function(x) exp(-.5*x)))
    colnames(thetasamps) <- c("Us","Ut")
    
  } else if (inherits(mod,'inla')) {
    # Get samples from INLA's marginals using the inverse CDF method
    betasamps <- vector(mode = 'list',length = length(mod$marginals.fixed))
    for (i in 1:length(mod$marginals.fixed)) betasamps[[i]] <- inla.qmarginal(runif(M),mod$marginals.fixed[[i]])
    betasamps <- Reduce(cbind,betasamps)
    colnames(betasamps) <- paste0('beta',1:ncol(betasamps))
    thetasamps <- vector(mode = 'list',length = length(mod$marginals.hyperpar))
    for (i in 1:length(mod$marginals.hyperpar)) thetasamps[[i]] <- inla.qmarginal(runif(M),mod$marginals.hyperpar[[i]])
    # Transform from PRECISION (INLA internally transforms from log-precision to precision) to STANDARD DEVIATION
    thetasamps <- Reduce(cbind,lapply(thetasamps,function(x) 1/sqrt(x)))
    colnames(thetasamps) <- c("Us","Ut")
  } else if (inherits(mod,'stanfit')) {
    modframe <- as.data.frame(mod)
    betasamps <- modframe[ ,grep('beta',colnames(modframe))]
    colnames(betasamps) <- paste0('beta',1:ncol(betasamps))
    # Transform from LOG-PRECISION to STANDARD DEVIATION
    thetasamps <- exp(-.5*modframe[ ,grep('theta',colnames(modframe))])
    colnames(thetasamps) <- c("Us","Ut")
  } else {
    stop(paste0("Unrecognized model type, ",class(mod),".\n"))
  }
  list(
    betasamps = betasamps,
    thetasamps = thetasamps
  )
}



# Obtain the data
cat("Reading in data, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

opioid_clean <- read_rds(datapath)

cat("Starting analysis, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")
set.seed(4564321)

cat("Fitting model, n = 1000000.\n")
models_n1000000 <- fit_all_models(sample_n(opioid_clean,1e06),domcmc = TRUE,tr=16)
cat("Fitting model, n = full.\n")
models_full <- fit_all_models(opioid_clean,domcmc = FALSE,tr=1)
cat("Done fitting models. Saving...\n")
save(models_n1000000,models_full,file = file.path(globalpath,"inla-elgm-full.RData"))
cat("Done saving models.\n")

## Create output table ----
# load and look
cat("Loading models...\n")
e <- new.env()
load(file.path(globalpath,"inla-elgm-full.RData"),envir = e)
cat("Done loading models.\n")
aghq_summary_n1000000 <- exp(summary(e$models_n1000000$aghq)$summarytable[ ,c(5,2,6)])
aghq_summary_nfull <- exp(summary(e$models_full$aghq)$summarytable[ ,c(5,2,6)])
inla_summary_n1000000 <- exp(e$models_n1000000$inla$internal.summary.hyperpar[ ,3:5])
inla_summary_nfull <- exp(e$models_full$inla$internal.summary.hyperpar[ ,3:5])

mcmc_n1000000 <- as.data.frame(e$models_n1000000$mcmc)[ ,c("theta[1]","theta[2]")]
mcmc_summary_n1000000 <- exp(data.frame(
  q025 = apply(mcmc_n1000000,2,quantile,probs = .025),
  q500 = apply(mcmc_n1000000,2,quantile,probs = .500),
  q975 = apply(mcmc_n1000000,2,quantile,probs = .975)
))
mcmc_summary_nfull <- data.frame(
  q025 = c(-1,-1),
  q500 = c(-1,-1),
  q975 = c(-1,-1)
)


# Compare (2.5,50,97.5) quantiles
quanttable <- data.frame(
  n = rep(rep(c(1e06,nrow(opioid_clean)),each=3),2),
  param = rep(paste0('theta',1:2),each=6),
  method = rep(c('aghq','inla','mcmc'),2),
  q025 = c(
    aghq_summary_n1000000[1,1],inla_summary_n1000000[1,1],mcmc_summary_n1000000[1,1],
    aghq_summary_nfull[1,1],inla_summary_nfull[1,1],mcmc_summary_nfull[1,1],
    aghq_summary_n1000000[2,1],inla_summary_n1000000[2,1],mcmc_summary_n1000000[2,1],
    aghq_summary_nfull[2,1],inla_summary_nfull[2,1],mcmc_summary_nfull[2,1]
  ),
  q500 = c(
    aghq_summary_n1000000[1,2],inla_summary_n1000000[1,2],mcmc_summary_n1000000[1,2],
    aghq_summary_nfull[1,2],inla_summary_nfull[1,2],mcmc_summary_nfull[1,2],
    aghq_summary_n1000000[2,2],inla_summary_n1000000[2,2],mcmc_summary_n1000000[2,2],
    aghq_summary_nfull[2,2],inla_summary_nfull[2,2],mcmc_summary_nfull[2,2]
  ),
  q975 = c(
    aghq_summary_n1000000[1,3],inla_summary_n1000000[1,3],mcmc_summary_n1000000[1,3],
    aghq_summary_nfull[1,3],inla_summary_nfull[1,3],mcmc_summary_nfull[1,3],
    aghq_summary_n1000000[2,3],inla_summary_n1000000[2,3],mcmc_summary_n1000000[2,3],
    aghq_summary_nfull[2,3],inla_summary_nfull[2,3],mcmc_summary_nfull[2,3]
  )
)

readr::write_csv(quanttable,file = file.path(globalpath,"inla-elgm-mcmc-quantile-table.csv"))
cat("Done.\n")
