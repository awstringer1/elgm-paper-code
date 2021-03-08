# Fit the disagg model using tmbstan
library(disaggregation)
library(tmbstan)
library(raster)

spatialdatapath <- "/home/alex/disaggregation_paper-master/"
savedatapath <- "/storage/tmbstan-model-20210205.RData"

## Load data ----

# This code copied from the software preprint https://arxiv.org/abs/2001.04847 accessed 2021/02/02
shapes <- shapefile(paste0(spatialdatapath,'data/shapes/mdg_shapes.shp')) 
population_raster <- raster(paste0(spatialdatapath,'data/population.tif'))
covariate_stack <- getCovariateRasters(paste0(spatialdatapath,'data/covariates'),
                                       shape = population_raster)

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
  ncores = 16
)

tau_u <- .1
tau_alpha <- .01
sigma_u <- 1
sigma_alpha <- .01
rho_u <- 3
rho_alpha <- .01

## BEGIN: taken directly from disaggregate::make_model_object()
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

limits <- sp::bbox(data$polygon_shapefile)
hypontenuse <- sqrt((limits[1,2] - limits[1,1])^2 + (limits[2,2] - limits[2,1])^2)
prior_rho <- hypontenuse/3

prior_sigma <- sd(data$polygon_data$response/mean(data$polygon_data$response))

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
## END: taken directly from disaggregate::make_model_object()

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


start <- Sys.time()
cat("Starting STAN, time =",format(start),"\n")
mcmc_out <- tmbstan(
  obj, 
  chains = 8, 
  iter = 10000, 
  warmup = 2000,
  cores = 8
)
end <- Sys.time() 
save(mcmc_out,file = savedatapath)
cat("Finished STAN, time:",format(difftime(end,start,units='secs')))

