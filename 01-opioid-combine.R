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
# Uncomment the next two lines if you don't want to use PARDISO when running this code:
inla.setOption(pardiso.license="/home/alex/sys/licenses/pardiso.lic")
stopifnot(inla.pardiso.check() == "SUCCESS: PARDISO IS INSTALLED AND WORKING")


# The number of iterations used in the MCMC time comparisons
# This is only for calculating the "number of equivalent iterations" below.
numitr <- 1000 

## Set paths ----

globalpath <- "/home/alex/data"
datapath <- file.path(globalpath,"opioid-admission-data-rds.rds")

aghqtimepath <- file.path(globalpath,"timeframe-20211022.csv")
inlatimepath <- file.path(globalpath,"inla-timeframe-20211022.csv")
mcmctimepath <- file.path(globalpath,"mcmc-timeframe-20211022.csv")

# Compile and load the TMB template. Source code included in the aghq package
tmbpath <- tempdir()
stopifnot(file.copy(system.file('extsrc/01_opioid.cpp',package='aghq'),tmbpath))
compile(file.path(tmbpath,"01_opioid.cpp"))
dyn.load(dynlib(file.path(tmbpath,"01_opioid")))

## END SETUP ##

## Fit all three methods to the same datasets ----
fit_all_models <- function(dat,domcmc=TRUE,tr=16) {
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
      int.strategy = 'ccd'
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
      iter = 1e05,
      warmup = 1e03,
      cores = 1 # The template is already parallelized, so don't do this.
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

# Obtain the data
cat("Reading in data, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")

opioid_clean <- read_rds(datapath)

cat("Starting analysis, time = ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n")
set.seed(4564321)
dat_1000 <- dplyr::sample_n(opioid_clean,1e03)
dat_10000 <- dplyr::sample_n(opioid_clean,1e04)
dat_100000 <- dplyr::sample_n(opioid_clean,1e05)
dat_1000000 <- dplyr::sample_n(opioid_clean,1e06)

cat("Fitting model, n = 1000.\n")
models_1000 <- fit_all_models(dat_1000)
cat("Fitting model, n = 10000.\n")
models_10000 <- fit_all_models(dat_10000)
cat("Fitting model, n = 100000.\n")
models_100000 <- fit_all_models(dat_100000)
cat("Fitting model, n = 1000000.\n")
models_1000000 <- fit_all_models(dat_1000000)
cat("Fitting model, n = full.\n")
models_full <- fit_all_models(opioid_clean,domcmc = FALSE,tr=1)

aghqenv <- new.env()
assign('model_obj_n1000',models_1000$aghq,envir = aghqenv)
assign('model_obj_n10000',models_10000$aghq,envir = aghqenv)
assign('model_obj_n100000',models_100000$aghq,envir = aghqenv)
assign('model_obj_n1000000',models_1000000$aghq,envir = aghqenv)
assign('model_obj_nfull',models_full$aghq,envir = aghqenv)

inlaenv <- new.env()
assign('model_obj_n1000',models_1000$inla,envir = inlaenv)
assign('model_obj_n10000',models_10000$inla,envir = inlaenv)
assign('model_obj_n100000',models_100000$inla,envir = inlaenv)
assign('model_obj_n1000000',models_1000000$inla,envir = inlaenv)
assign('model_obj_nfull',models_full$inla,envir = inlaenv)

mcmcenv <- new.env()
assign('model_obj_n1000',models_1000$mcmc,envir = mcmcenv)
assign('model_obj_n10000',models_10000$mcmc,envir = mcmcenv)
assign('model_obj_n100000',models_100000$mcmc,envir = mcmcenv)
assign('model_obj_n1000000',models_1000000$mcmc,envir = mcmcenv)
assign('model_obj_nfull',models_full$mcmc,envir = mcmcenv)

cat("Finished fitting models. Moving on to summaries.\n")

## Compute summaries ----

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

# Now actually obtain the summaries
summarylist <- list()
for (nm in c('aghq','inla','mcmc')) {
  cat(paste0('nm = ',nm,":"),sep = '')
  summarylist[[nm]] <- list()
  for (nn in names(aghqenv)) {
    cat(paste0(' nn = ',nn,","),sep='')
    ee <- get(paste0(nm,'env'),envir = .GlobalEnv)
    tmp <- tryCatch(get(nn,envir = ee),error = function(e) 'next')
    if (!inherits(tmp,'condition')) if (!inherits(tmp,'list')) summarylist[[nm]][[nn]] <- get_samples(tmp)
  }
  cat(".\n")
}

## Point estimates ----
# Not reported in paper, used as a quick check.
coeflists <- lapply(summarylist,function(lst) lapply(lst,function(itm) Reduce(c,lapply(itm,function(x) apply(x,2,mean)))))
coefframe <- dplyr::bind_rows(lapply(coeflists,dplyr::bind_rows,.id = "model"),.id = "method") %>% mutate(model = stringr::str_remove(model,'model_obj_'))

## KS ----
# Compute the KS between each of INLA and AGHQ with MCMC, and also each other, for each sample size
# This is Table 2
nms <- names(summarylist[[1]])
nm3 <- names(summarylist[[3]])
kstable <- list()
for (nm in nms) {
  aghqframe <- Reduce(cbind,summarylist$aghq[[nm]])
  inlaframe <- Reduce(cbind,summarylist$inla[[nm]])
  if (nm %in% nm3) mcmcframe <- Reduce(cbind,summarylist$mcmc[[nm]])
  aghqks <- inlaks <- aghqinlaks <- numeric(ncol(aghqframe))
  for (i in 1:length(aghqks)) {
    # Suppress warnings about approximate p-values
    if (nm %in% nm3) {
      suppressWarnings(aghqks[i] <- ks.test(aghqframe[ ,i],mcmcframe[ ,i])$statistic)
      suppressWarnings(inlaks[i] <- ks.test(inlaframe[ ,i],mcmcframe[ ,i])$statistic)
    } else{
      aghqks[i] <- inlaks[i] <- -1
    }
    suppressWarnings(aghqinlaks[i] <- ks.test(inlaframe[ ,i],aghqframe[ ,i])$statistic)
  }
  kstable <- c(kstable,list(data.frame(n = nm,var = c(paste0('beta',1:8),'Us','Ut'),aghq = aghqks,inla = inlaks,aghqinla = aghqinlaks)))
}
kstable <- Reduce(bind_rows,kstable) %>% mutate(n = stringr::str_remove(n,'model_obj_'))

readr::write_csv(kstable,file.path(globalpath,"ks-table-20211101.csv"))

## Time Comparisons ----
# Read in the tables of timings
aghqtime <- readr::read_csv(aghqtimepath,col_names = TRUE,col_types = c('nnn')) %>% dplyr::mutate(method = 'aghq')
inlatime <- readr::read_csv(inlatimepath,col_names = TRUE,col_types = c('nnn')) %>% dplyr::mutate(method = 'inla')
mcmctime <- readr::read_csv(mcmctimepath,col_names = TRUE,col_types = c('nnn')) %>% dplyr::mutate(method = 'mcmc')

# Compute mean and standard deviation
timecompare <- dplyr::bind_rows(aghqtime,inlatime,mcmctime) %>%
  group_by(method,n) %>%
  summarize(meantime = mean(time),sdtime = sd(time))

# Compute number of equivalent iterations
# numitr defined at top.
mcmctimes <- dplyr::filter(timecompare,method=='mcmc')
equiviter <- dplyr::inner_join(filter(timecompare,method!='mcmc'),mcmctimes,by='n') %>%
  dplyr::mutate(equiviter = (numitr / meantime.y) * meantime.x) %>%
  dplyr::select(method=method.x,n,equiviter)

# This is Table 1
timecompare <- dplyr::left_join(timecompare,equiviter,by=c('method','n'))

readr::write_csv(timecompare,file.path(globalpath,"timecompare-table-20211101.csv"))

cat("Done.\n")