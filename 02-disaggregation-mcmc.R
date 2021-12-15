### Spatial Disaggregation Example ###
# This code runs MCMC using tmbstan for example 6.1
# The data are downloaded and everything is saved to temp storage
# You should not have to change anything in order to source this script in a clean
# session.

## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'disaggregation',
  'raster',
  'tmbstan'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    if (pkg == 'disaggregation') {
      if (!require('devtools',character.only = TRUE,quietly = TRUE)) install.packages('devtools')
      devtools::install_github('aknandi/disaggregation')
      require(pkg,character.only = TRUE,quietly = TRUE)
    } else {
      install.packages(pkg)
      require(pkg,character.only = TRUE,quietly = TRUE)
    }
  }
}

## Global variables

# Path to spatial data from disaggregation package
# https://github.com/aknandi/disaggregation_paper

spatialdatapath <- "https://github.com/aknandi/disaggregation_paper/blob/master/data/"
populationdatapath <- paste0(spatialdatapath,"population.tif?raw=true")
shpfl <- c("dbf","prj","shp","shx")
shapefiledatapath <- paste0(spatialdatapath,"shapes/mdg_shapes.",shpfl,"?raw=TRUE")
covfl <- c("EVI","Elevation","LSTmean")
covdatapath <- paste0(spatialdatapath,"covariates/",covfl,".tif?raw=TRUE")

# Path to save mcmc results
# Creating the mcmc results took ~29 hours on my hardware.
mcmcdir <- file.path(tempdir(),"mcmc")
if (!dir.exists(mcmcdir)) dir.create(mcmcdir)
## Load data ----

# Download the data into temp storage
datadir <- file.path(tempdir(),"data")
if (!dir.exists(datadir)) dir.create(datadir)
populationdatafile <- file.path(datadir,"populationdata.tif")
shapedir <- file.path(datadir,"shapefiles")
if (!dir.exists(shapedir)) dir.create(shapedir)
shapefiledatafile <- c(
  file.path(shapedir,"mdg_shapes.dbf"),
  file.path(shapedir,"mdg_shapes.prj"),
  file.path(shapedir,"mdg_shapes.shp"),
  file.path(shapedir,"mdg_shapes.shx")
)
covariatedir <- file.path(datadir,"covariates")
if (!dir.exists(covariatedir)) dir.create(covariatedir)
covdatafile <- c(
  file.path(covariatedir,"EVI.tif"),
  file.path(covariatedir,"Elevation.tif"),
  file.path(covariatedir,"LSTmean.tif")
)

# Download
# Population data
tryCatch(
  download.file(populationdatapath,populationdatafile),
  warning = function(w) cat(paste0("Could not download population data from ",populationdatapath," to ",populationdatafile,". Received the following warning message: ",w,".\n")),
  warning = function(e) cat(paste0("Could not download population data from ",populationdatapath," to ",populationdatafile,". Received the following error message: ",e,".\n"))
)
stopifnot(file.exists(populationdatafile))
# Shape files
for (i in 1:length(shpfl)) {
  tmp1 <- shapefiledatapath[i]
  tmp2 <- shapefiledatafile[i]
  tryCatch(
    download.file(tmp1,tmp2),
    warning = function(w) cat(paste0("Could not download shapefile from ",tmp1," to ",tmp2,". Received the following warning message: ",w,".\n")),
    warning = function(e) cat(paste0("Could not download shapefile from ",tmp1," to ",tmp2,". Received the following error message: ",e,".\n"))
  )
  stopifnot(file.exists(tmp2))
}
# Covariates
for (i in 1:length(covfl)) {
  tmp1 <- covdatapath[i]
  tmp2 <- covdatafile[i]
  tryCatch(
    download.file(tmp1,tmp2),
    warning = function(w) cat(paste0("Could not download shapefile from ",tmp1," to ",tmp2,". Received the following warning message: ",w,".\n")),
    warning = function(e) cat(paste0("Could not download shapefile from ",tmp1," to ",tmp2,". Received the following error message: ",e,".\n"))
  )
  stopifnot(file.exists(tmp2))
}

# Final check with informative error message if data aren't found in the correct location
# needed for reading in, in the following step
for (i in 1:length(shpfl)) {
  if (!file.exists(shapefiledatafile[i])) stop(paste0("Could not find shapefile ",shapefiledatafile[i]," which is needed. There should have been an error in the code above that downloads the data, please check.\n"))
}
for (i in 1:length(covfl)) {
  if (!file.exists(covdatafile[i])) stop(paste0("Could not find covariate raster ",covdatafile[i]," which is needed. There should have been an error in the code above that downloads the data, please check.\n"))
}
if (!file.exists(populationdatafile)) stop(paste0("Could not find population data raster ",populationdatafile," which is needed. There should have been an error in the code above that downloads the data, please check.\n"))


# This code copied from the software preprint https://arxiv.org/abs/2001.04847 accessed 2021/02/02
shapes <- shapefile(shapefiledatafile[3]) 
population_raster <- raster(populationdatafile)
covariate_stack <- getCovariateRasters(covariatedir,shape = population_raster)

# Takes a minute or two
dis_data <- prepare_data(
  polygon_shapefile = shapes, 
  covariate_rasters = covariate_stack,
  aggregation_raster = population_raster,
  mesh.args = list(max.edge = c(0.7, 8),
                   cut = 0.05,
                   offset = c(1, 2)),
  id_var = 'ID_2', 
  response_var = 'inc', 
  na.action = TRUE, 
  ncores = 8
)

## Fit their model ----

# Prior params, for later too

tau_u <- .1
tau_alpha <- .01
sigma_u <- 1
sigma_alpha <- .01
rho_u <- 3
rho_alpha <- .01

# This is the package code to fit the model, but it does not include
# beta in the marginal laplace, so is not directly usable.
# system.time(
# fitted_model <- disag_model(
#   data = dis_data, 
#   iterations = 1000,
#   family = 'poisson',
#   link = 'log',
#   priors = list(
#     priormean_intercept = 0,
#     priorsd_intercept = 2, 
#     priormean_slope = 0.0, 
#     priorsd_slope = 0.4, 
#     prior_rho_min = rho_u, 
#     prior_rho_prob = rho_alpha, 
#     prior_sigma_max = sigma_u, 
#     prior_sigma_prob = sigma_alpha)
#   )
# )
# Instead, replicate their template creation and then change the 'random' option.


## Prepare model object ----
# Need the TMB function object
# Will fit their model using their code, but for ELGM we need the regression
# coefficients to be part of the marginal laplace approx
# Here I make these minor modifications to their make_model_object() function.

## BEGIN: taken directly from disaggregate::make_model_object() ##
data = dis_data
family = 'poisson'
link = 'log'
priors = list( # Combination of priors they used plus default priors in the source code
  priormean_intercept = 0,
  priorsd_intercept = 2, 
  priormean_slope = 0.0, 
  priorsd_slope = 0.4, 
  prior_rho_min = rho_u, 
  prior_rho_prob = rho_alpha, 
  prior_sigma_max = sigma_u, 
  prior_sigma_prob = sigma_alpha,
  prior_iideffect_sd_max = tau_u,
  prior_iideffect_sd_prob = tau_alpha)

field <- iid <- TRUE

family_id <- 2 # Poisson
link_id <- 1 # Log link

nu = 1
spde <- (INLA::inla.spde2.matern(data$mesh, alpha = nu + 1)$param.inla)[c("M0", "M1", "M2")]	
Apix <- INLA::inla.mesh.project(data$mesh, loc = data$coordsForFit)$A
n_s <- nrow(spde$M0)

cov_matrix <- as.matrix(data$covariate_data[, -c(1:2)])
cov_matrix <- t(apply(cov_matrix, 1,as.numeric))

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)),
                   log_tau_gaussian = 8,
                   iideffect = rep(0, nrow(data$polygon_data)),
                   iideffect_log_tau = 1,
                   log_sigma = 0,
                   log_rho = 4,
                   nodemean = rep(0, n_s))

input_data <- list(x = cov_matrix,
                   aggregation_values = data$aggregation_pixels,
                   Apixel = Apix,
                   spde = spde,
                   startendindex = data$startendindex,
                   polygon_response_data = data$polygon_data$response,
                   response_sample_size = data$polygon_data$N,
                   family = family_id,
                   link = link_id,
                   nu = nu,
                   field = as.integer(field),
                   iid = as.integer(iid))

input_data <- c(input_data, priors)

tmb_map <- list(log_tau_gaussian = as.factor(NA))
## END: taken directly from disaggregate::make_model_object() ##

# CHANGE: include regression coefficients in the Laplace approx
random_effects <- c('nodemean','iideffect','intercept','slope')

obj <- TMB::MakeADFun(
  data = input_data, 
  parameters = parameters,
  map = tmb_map,
  random = random_effects,
  hessian = TRUE,
  silent = TRUE, # Always do silent = TRUE
  DLL = "disaggregation")

## Start fit ----

start <- Sys.time()
cat("Starting STAN, time =",format(start),"\n")
mcmc_out <- tmbstan(
  obj, 
  chains = 16, 
  iter = 5000, 
  warmup = 1000,
  cores = 16
)
end <- Sys.time() 
save(mcmc_out,file = file.path(mcmcdir,"mcmc-results-20211022.RData"))
cat("Finished STAN, time:",format(difftime(end,start,units='secs')))
cat("Done.\n")
