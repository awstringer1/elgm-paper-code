### Cox PH example ###
# Reproduce example 6.2
# This script should be sourcable in a clean session without modification.
# NOTE: UK border data for plots are downloaded using raster::get_data
# Sometimes the server is down. In this case, the script fits the model but doesn't
# compute any summaries.

## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'Matrix',
  'aghq',
  'sp',
  'geostatsp',
  'raster',
  'INLA' # For Leuk data only
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    if (pkg == "INLA") {
      install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
      # Don't load it, only need for Leuk data
    } else {
      install.packages(pkg)
      require(pkg,character.only = TRUE,quietly = TRUE)
    }
  }
}



## Global Variables ----
options(mc.cores = 16)
globalpath <- tempdir()
datapath <- file.path(globalpath,"data")
spatialdatapath <- file.path(datapath,"aggspatial/")
figurepath <- file.path(globalpath,"figures")
if (!dir.exists(datapath)) dir.create(datapath)
if (!dir.exists(spatialdatapath)) dir.create(spatialdatapath)
if (!dir.exists(figurepath)) dir.create(figurepath)

datestamp <- "20210211"
version <- "v1"

savepath <- file.path(datapath,paste0("coxph-model-",datestamp,"-",version,".RData"))
savebrickpath <- file.path(datapath,paste0("coxph-brick-",datestamp,"-",version,".RData"))

domodel <- TRUE # TRUE: fit model, FALSE: load from disk
dobrick <- TRUE # TRUE: simulate fields, FALSE: load from disk
dosummaries <- TRUE # FALSE: don't bother with computing summaries.

## Download the UK border data ----
# Sometimes this fails, because the server is down.
# In that case, do not do the summaries
ukBorderLL <- tryCatch(raster::getData("GADM", country='GBR', level=3),error = function(e) e)
if (inherits(ukBorderLL,'condition')) {
  cat(paste0("Could not download the UK border data. Setting dobrick and dosummaries to FALSE.\n"))
  dobrick <- dosummaries <- FALSE
}

PLOTTEXTSIZE <- 28

set.seed(573892)

# Prior distribution parameters
sigma_u <- 1
sigma_alpha <- .5
rho_u <- 20 * 1000
rho_alpha <- .5
beta_prec <- .001

## Load data ----
data(Leuk,package = "INLA")
leuk <- as_tibble(Leuk)

# Centre variables, add a small amount of normal noise for removing ties,
# sort data in increasing order of survival with censored obs at the end

agecentre <- mean(leuk$age)
wbccentre <- mean(leuk$wbc)
tpibin <- round(leuk$tpi,1)
tpicentre <- 0

dataformodelling <- leuk %>%
  mutate(time = time + rnorm(length(leuk$time),sd = .00001),
         age = age - agecentre,
         wbc = wbc - wbccentre,
         tpi = tpibin - tpicentre,
  ) %>%
  arrange(desc(cens),time)

# Index of centre value
tpiidx <- which(sort(unique(dataformodelling$tpi)) == 0)

# Add the spatial component.
# Make a spatial points dataframe containing the point locations of each incidence
# Then overlay a raster which covers northeast england
PROJTOUSE <- mapmisc::omerc(c(-3.055,53.365),angle=0)
crstouse <- CRS("+init=epsg:27700")

pointsdata <- SpatialPointsDataFrame(
  coords = 90*1000*dplyr::select(dataformodelling,xcoord,ycoord),
  data = dplyr::select(dataformodelling,-xcoord,-ycoord),
  proj4string = PROJTOUSE
)
pointsdata <- spTransform(pointsdata,crstouse)
      
## END SETUP ##

## Model setup ----

# Design matrix
Amat <- Diagonal(n = nrow(dataformodelling),x = 1) 
Xmat <- sparse.model.matrix(time ~ age + sex + wbc + tpi,data = dataformodelling)

# Censoring indicator- 1 means NOT censored
# Remove the first observation (guaranteed not to be censored by sorting), as likelihood
# depends only on differences from this obs
censor <- dataformodelling$cens[-1] # 1 == not censored, confusing but more useful.

# Differencing
create_diff_matrix <- function(n) {
  cbind(Matrix(1,n-1,1),Diagonal(n-1,-1))
}

# Differenced design matrix

diffmat <- create_diff_matrix(nrow(dataformodelling))
Zmat <- diffmat %*% cbind(Amat,Xmat) 

make_delta <- function(W) {
  as.numeric(Zmat %*% cbind(W))
}

# Dimensions
p <- ncol(Xmat)
d <- ncol(Amat)
Wd <- ncol(Zmat)
stopifnot(p+d == Wd)

# Likelihood and derivatives

compute_one_denominator <- function(delta,i) {
  # All of the likelihood quantities require that denominator
  # vector for each observation. It's a cumulative sum. Write
  # one function that computes it efficiently.
  # delta: vector of length n
  # i: index of denominator you want
  n <- length(delta)
  dd <- delta[i] - delta[i:n]
  exp(matrixStats::logSumExp(dd)) - 1
}

compute_denominator <- function(delta) {
  map(1:length(delta),~compute_one_denominator(delta,.x)) %>% reduce(c)
}

log_likelihood <- function(W) {
  delta <- make_delta(W)
  denom <- compute_denominator(delta)
  -sum(censor * log(1 + denom))
}

grad_log_likelihood_one <- function(W,i) {
  delta <- make_delta(W)
  n <- length(delta)
  if (censor[i] == 0) return(sparseVector(0,0,n))
  denom <- compute_one_denominator(delta,i)
  dd <- delta[i] - delta[i:n]
  out <- exp(dd) / (1 + denom)
  out <- c(rep(0,i-1),out)
  out[i] <- out[i] - 1
  out
}

grad_log_likelihood_subset <- function(W,I) {
  # I: subset of 1...n
  map(I,~grad_log_likelihood_one(W,.x)) %>% reduce(~.x + .y)
}

grad_log_likelihood <- function(W) grad_log_likelihood_subset(W,which(censor==1))


make_hess_vec <- function(delta,i) {
  # Make the vector that is used to create the hessian
  n <- length(delta)
  denom <- compute_one_denominator(delta,i)
  if (censor[i] == 0) return(sparseVector(0,0,n))
  dd <- delta[i] - delta[i:n]
  out <- exp(dd) / (1 + denom)
  c(rep(0,i-1),out)
}


hessian_log_likelihood <- function(W) {
  delta <- make_delta(W)
  n <- length(delta)
  gg <- map(which(censor==1),~make_hess_vec(delta,.x))
  diag(as.numeric(reduce(gg,~.x+.y))) - tcrossprod(gg %>% reduce(cbind))
}

# Prior

Q_matrix <- function(theta) {
  # theta = log(sigma), log(rho)
  theta <- as.numeric(unname(theta))
  # Matern
  mm <- geostatsp::matern(
    pointsdata,
    param = c("variance" = exp(2 * theta[1]),"range" = exp(theta[2]),"shape" = 1),
    type = "precision"
  )
  
  bb <- beta_prec * diag(p)
  
  rbind(
    cbind(mm,Matrix(0,nrow = nrow(mm),ncol = p,sparse = FALSE)),
    cbind(Matrix(0,nrow = p,ncol = ncol(mm),sparse = FALSE),bb)
  )
}

logsigmalogprior <- function(theta,sigma0,prior_alpha) {
  # theta = log(sigma)
  alpha2 <- prior_alpha
  lambda2 <- -log(alpha2) / sigma0
  log(lambda2) - exp(theta) * lambda2 + theta
}

logrhologprior <- function(theta,rho0,prior_alpha) {
  # theta = log(rho)
  d <- 2
  alpha1 <- prior_alpha
  lambda1 <- -log(alpha1) * rho0^(d/2)
  log(d/2) + log(lambda1) - (d/2)*theta - lambda1 * exp(-theta * d/2)
}

logprior_theta <- function(theta) 
  logsigmalogprior(theta[1],sigma0 = sigma_u,prior_alpha = sigma_alpha) + 
  logrhologprior(theta[2],rho0 = rho_u,prior_alpha = rho_alpha)

log_prior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  -(1/2) * as.numeric(crossprod(W,crossprod(Q,W))) +(1/2)*as.numeric(determinant(Q,logarithm = TRUE)$modulus)
}

log_posterior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  log_prior_W(W,theta,Q) + log_likelihood(W) + logprior_theta(theta)
}

grad_log_posterior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  as.numeric(-Q %*% cbind(W) + t(Zmat) %*% grad_log_likelihood(W))
}

H_matrix <- function(W,theta) {
  C <- hessian_log_likelihood(W)
  Q <- Q_matrix(theta)
  
  Q + crossprod(Zmat,crossprod(C,Zmat))
}

## Fit model ----

ff <- list(
  fn = log_posterior_W,
  gr = grad_log_posterior_W,
  he = function(W,theta) -1 * H_matrix(W,theta)
)

if (domodel) {
  tm <- Sys.time()
  cat("Fitting model, time = ",format(tm),"\n")
  coxphmod <- marginal_laplace(
    ff,
    k = 3,
    # theta starting values from Brown (2015), geostatsp/diseasemapping software paper
    startingvalue = list(W = rep(0,Wd),theta = c(-1,log(50*1000))),
    control = default_control_marglaplace(method = 'BFGS',inner_method = 'trust')
  )
  cat("Finished model, took ",format(difftime(Sys.time(),tm,units = 'secs')),"\n")
  save(coxphmod,file = savepath)
} else {
  load(savepath)
}

## Sample from the posterior ----

set.seed(87964382)
posterior_samples <- sample_marginal(coxphmod,1e03)

if (dosummaries) {
  cat("Posterior summaries, time = ",format(tm),"\n")
  # Parameter table
  # Betas
  betaidx <- (Wd-p+2):Wd # Note: first index is intercept, not actually estimable here
  betasamps <- posterior_samples$samps[betaidx, ]
  
  betamean <- apply(betasamps,1,mean)
  betasd <- apply(betasamps,1,sd)
  beta2.5 <- apply(betasamps,1,quantile,probs = .025)
  beta97.5 <- apply(betasamps,1,quantile,probs = .975)
  
  # Sigmas
  sigma_and_rho_means <- compute_moment(coxphmod$normalized_posterior,exp)
  sigma_and_rho_sd <- sqrt(compute_moment(coxphmod$normalized_posterior,function(x) ((exp(x) - sigma_and_rho_means)^2)))
  sigma_quants <- exp(compute_quantiles(coxphmod$marginals[[1]]))
  rho_quants <- exp(compute_quantiles(coxphmod$marginals[[2]]))
  
  coeftable <- tibble(
    mean = c(betamean,sigma_and_rho_means),
    sd = c(betasd,sigma_and_rho_sd),
    q2.5 = c(beta2.5,sigma_quants[1],rho_quants[1]),
    q97.5 = c(beta97.5,sigma_quants[2],rho_quants[2])
  )
  
  readr::write_csv(coeftable,file.path(figurepath,paste0("coxphcoeftable-",datestamp,"-",version,".csv")))
  
  knitr::kable(
    coeftable,
    digits = 3,
    format = 'latex'
  )
  
  
  # Plot the marginals with priors
  sigmaprior <- function(sigma,u,alpha) dexp(sigma,-log(alpha)/u) # Same for tau as well
  rhoprior <- function(rho,u,alpha) (1/rho^2) * dexp(1/rho,-log(alpha) * u)
  
  sigmapdf <- compute_pdf_and_cdf(coxphmod$marginals[[1]],list(totheta = log,fromtheta = exp),seq(-10,-.2,by=.01))
  rhopdf <- compute_pdf_and_cdf(coxphmod$marginals[[2]],list(totheta = log,fromtheta = exp),seq(0,12.5,by=.01))
  
  
  # plot
  PLOTTEXTSIZE <- 28
  
  sigmaplot <- sigmapdf %>%
    ggplot(aes(x = transparam,y = pdf_transparam)) +
    theme_classic() +
    geom_line() +
    stat_function(fun = sigmaprior,args = list(u = sigma_u,alpha = sigma_alpha),linetype = 'dashed') +
    # labs(x = expression(sigma),y = "Density") +
    labs(x = "",y = "") +
    theme(text = element_text(size = PLOTTEXTSIZE))
  
  rhoplot <- rhopdf %>%
    ggplot(aes(x = transparam,y = 1e05*pdf_transparam)) +
    theme_classic() +
    geom_line() +
    geom_line(data = tibble(x = seq(0.1,2e05,length.out = 1e05),y = 1e05*rhoprior(x,rho_u,rho_alpha)),aes(x = x,y = y),linetype = "dashed") +
    scale_x_continuous(breaks = seq(0,1e06,by = 2.5e04),labels = function(x) x/1000) +
    theme(text = element_text(size = PLOTTEXTSIZE)) +
    coord_cartesian(xlim = c(0,1.5e05)) +
    labs(x = "",y = "")
  
  WIDTH <- HEIGHT <- 7
  ggsave(
    file.path(figurepath,paste0("coxphsigmapostplot-",datestamp,"-",version,".pdf")),
    sigmaplot,
    width = WIDTH,height = HEIGHT
  )
  ggsave(
    file.path(figurepath,paste0("coxphrhopostplot-",datestamp,"-",version,".pdf")),
    rhoplot,
    width = WIDTH,height = HEIGHT
  )
  
}

# Simulate the spatial fields
simulate_spatial_fields <- function(U,
                                    theta,
                                    pointsdata,
                                    resolution = list(nrow = 100,ncol = 100)) {
  # U: matrix of samples, each column is a sample
  # theta: tibble of theta values
  # Draw from U*|U
  fieldlist <- vector(mode = 'list',length = nrow(theta))
  for (i in 1:length(fieldlist)) {
    fielddat <- pointsdata
    fielddat@data <- data.frame(w = as.numeric(U[ ,i]))
    
    # Back-transform the Matern params
    sig <- exp(theta$theta1[i])
    rho <- exp(theta$theta2[i])
    # Simulate from the two fields
    capture.output({
      fieldlist[[i]] <- geostatsp::RFsimulate(
        model = c("variance" = sig^2,"range" = rho,"shape" = 2),
        data = fielddat,
        x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol),
        n = 1
      )
    })
  }
  brick(fieldlist)
}

if (dobrick) {
  tm <- Sys.time()
  cat("Doing brick, time = ",format(tm),"\n")
  # Do it for only 100, because it takes a lot of time
  set.seed(8796)
  simstodo <- sample(ncol(posterior_samples$samps),100,replace = FALSE)
  fieldsbrick <- simulate_spatial_fields(
    U = posterior_samples$samps[1:d,simstodo],
    theta = posterior_samples$theta[simstodo, ],
    pointsdata = pointsdata,
    resolution = list(nrow = 400,ncol = 200)
  )
  cat("Finished brick, took ",format(difftime(Sys.time(),tm,units = 'secs')),"\n")
  save(fieldsbrick,file = savebrickpath)
} else {
  load(savebrickpath)
}

## Plots ----

if (dosummaries) {
  ukBorder = spTransform(ukBorderLL, projection(pointsdata))
  pointsinpoly <- pointsdata %over% ukBorder
  pointsinpolyID <- unique(pointsinpoly$GID_2)
  ukBorder <- ukBorder[ukBorder$GID_2 %in% pointsinpolyID, ]
  # Get the outer border
  ukBorderouter <- rgeos::gUnaryUnion(ukBorder)
  
  
  simfieldsmean <- mean(exp(fieldsbrick))
  simfieldsexceedence <- mean(fieldsbrick > log(1.2))
  
  # MEAN
  plotraster <- simfieldsmean
  
  predcols <- mapmisc::colourScale(
    plotraster,
    breaks = quantile(values(plotraster),probs = (0:9)/9),
    style = "fixed",
    col = "Spectral",
    rev = TRUE,
    dec = -log10(0.05)
  )
  
  colvals <- 100
  bigvalues <- quantile(values(plotraster),probs = (0:(colvals-1))/(colvals-1))
  plotraster <- mask(plotraster, ukBorder)
  
  
  pdf(file.path(figurepath,paste0("coxphspatialplot-",datestamp,"-",version,".pdf")),width = 7,height = 7)
  mapmisc::map.new(pointsdata)
  plot(plotraster,
       col = predcols$col,
       breaks = predcols$breaks,
       legend=FALSE, add=TRUE)
  plot(plotraster, col=predcols$col, breaks=predcols$breaks, legend=FALSE, add=TRUE)
  plot(ukBorder, add=TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
  plot(ukBorderouter,add = TRUE)
  points(pointsdata,pch = ".")
  mapmisc::legendBreaks('topright', predcols, cex=1.5, bty='n',inset=0)
  dev.off()
  
  # EXCEEDENCE PROBABILITIES
  
  plotraster <- simfieldsexceedence
  
  predcols <- mapmisc::colourScale(
    plotraster,
    breaks = c(0,0.05,0.1,0.15,0.2,0.3,0.5,0.65),
    style = "fixed",
    col = "Spectral",
    rev = TRUE
  )
  
  plotraster <- mask(plotraster, ukBorder)
  
  pdf(file.path(figurepath,paste0("coxphexceedanceplot-",datestamp,"-",version,".pdf")),width = 7,height = 7)
  mapmisc::map.new(pointsdata)
  plot(plotraster,
       col = predcols$col,
       breaks = predcols$breaks,
       legend=FALSE, add=TRUE)
  plot(ukBorder, add=TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
  plot(ukBorderouter,add = TRUE)
  points(pointsdata,pch = ".")
  mapmisc::legendBreaks('topright', predcols, cex=1.5, bty='n')
  dev.off()
}
cat(paste0("Done. You can go to ",figurepath," to see the output. Thanks!\n"))