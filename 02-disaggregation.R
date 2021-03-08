### Spatial Disaggregation Example ###

## Load libraries ----
library(aghq)
library(disaggregation)
library(raster)
library(tidyverse)

## Global variables

# Path to spatial data from disaggregation package
# https://github.com/aknandi/disaggregation_paper
spatialdatapath <- "/storage/disaggregation_paper-master/"

# Path to mcmc results
mcmcpath <- "/storage/tmbstan-model-20210205.RData"

# Path to save figures
figurepath <- "/storage/phd/projects/no-epsilon/jcgs/paper/figures/"
datestamp <- "20210210"
version <- "v1"


PLOTTEXTSIZE <- 28
NLPLOTTEXTSIZE <- 18


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
# beta in the marginal laplace, so is not directly comparable to our approach
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


## Prepare model object ----
# Need the TMB function object
# Will fit their model using their code, but for ELGM we need the regression
# coefficients to be part of the marginal laplace approx
# Further, it's better here if we don't use TMB's Laplace, since we have our own
# Here I make these minor modifications to their make_model_object() function.

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

## Fit the model using a Laplace approx to all params ----
logpost_fn <- function(theta) -1 * obj$fn(theta)
logpost_gr <- function(theta) -1 * obj$gr(theta)
logpost_he <- function(theta) -1 * numDeriv::jacobian(obj$gr,theta)
ff <- list(fn = logpost_fn,gr = logpost_gr,he = logpost_he)

system.time({
  opt <- aghq::optimize_theta(ff,startingvalue = rep(0,3),control = list(method = 'trust'))
  sdr <- TMB::sdreport(obj)
  sdsummary <- summary(sdr,select = 'fixed')
}) # 95 seconds


## Fit the model using AGHQ ----
# theta = log(tau), log(sigma), log(rho)

system.time(
  aghqmodel <- aghq::aghq(ff,k=7,startingvalue = opt$mode,optresults = opt)
) # k = 7, pre-optimized: 1642 seconds

## Inference for theta ----
summary(aghqmodel)
# Transformations
taupdf <- compute_pdf_and_cdf(aghqmodel$marginals[[1]],
                              transformation = list(
                                totheta = function(x) -2*log(x),
                                fromtheta = function(x) exp(-.5*x)))
sigmapdf <- compute_pdf_and_cdf(aghqmodel$marginals[[2]],transformation = list(totheta = log,fromtheta = exp))
rhopdf <- compute_pdf_and_cdf(aghqmodel$marginals[[3]],transformation = list(totheta = log,fromtheta = exp))

# Priors
sigmaprior <- function(sigma,u,alpha) dexp(sigma,-log(alpha)/u) # Same for tau as well
rhoprior <- function(rho,u,alpha) (1/rho^2) * dexp(1/rho,-log(alpha) * u)

# plot

tauplot <- taupdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_bw() +
  geom_line() +
  stat_function(fun = sigmaprior,args = list(u = tau_u,alpha = tau_alpha),linetype = 'dashed') +
  stat_function(fun = dlnorm,args = list(meanlog = sdsummary['iideffect_log_tau','Estimate'],sdlog = sdsummary['iideffect_log_tau','Std. Error']),linetype = 'dotdash') +
  labs(x = expression(tau),y = "Density") +
  theme(text = element_text(size = PLOTTEXTSIZE))

sigmaplot <- sigmapdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_bw() +
  geom_line() +
  stat_function(fun = sigmaprior,args = list(u = sigma_u,alpha = sigma_alpha),linetype = 'dashed') +
  stat_function(fun = dlnorm,args = list(meanlog = sdsummary['log_sigma','Estimate'],sdlog = sdsummary['log_sigma','Std. Error']),linetype = 'dotdash') +
  labs(x = expression(sigma),y = "Density") +
  theme(text = element_text(size = PLOTTEXTSIZE))

rhoplot <- rhopdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_bw() +
  geom_line() +
  stat_function(fun = rhoprior,args = list(u = rho_u,alpha = rho_alpha),linetype = 'dashed') +
  stat_function(fun = dlnorm,args = list(meanlog = sdsummary['log_rho','Estimate'],sdlog = sdsummary['log_rho','Std. Error']),linetype = 'dotdash') +
  labs(x = expression(rho),y = "Density") +
  theme(text = element_text(size = PLOTTEXTSIZE)) +
  coord_cartesian(xlim = c(0,5),ylim = c(0,1))


## Inference for W ----
# Get a list of modes and Hessians at the quadrature points

distinctthetas <- as_tibble(
  aghqmodel$normalized_posterior$nodesandweights
) %>%
  select(contains('theta'))
modesandhessians <- distinctthetas %>% 
  tibble::add_column(mode = vector(mode = "list", length = nrow(distinctthetas)), 
                     H = vector(mode = "list", length = nrow(distinctthetas)))

system.time({
  for (i in 1:nrow(distinctthetas)) {
    cat('Iteration',i,'of',nrow(distinctthetas),'\n')
    # Get the theta
    theta <- as.numeric(modesandhessians[i,c('theta1','theta2','theta3')])
    # Set up the mode and hessian of the random effects. This happens when you run
    # the TMB objective with a particular theta
    obj$fn(theta)
    # Now pull the mode and hessian
    modesandhessians[i,'mode'] <- list(list(obj$env$last.par[!(names(obj$env$last.par) %in% c("iideffect_log_tau","log_sigma","log_rho"))]))
    modesandhessians[i,'H'] <- list(list(obj$env$spHess(obj$env$last.par,random = TRUE)))
  }
})

aghqmodel$modesandhessians <- modesandhessians
samps <- sample_marginal(aghqmodel,1e02)

## Plot the maps ----

# Get the background etc
MGBorderLL = raster::getData("GADM", country='MG', level=3) # Regions
MGBorder = spTransform(MGBorderLL, projection(dis_data$polygon_shapefile))
# Get the outer border
MGBorderouter <- rgeos::gUnaryUnion(MGBorder)

# Plot function
plotmap <- function(rast,brk,saveit="",inset = 0) {
  rast <- raster::mask(rast,MGBorder)
  predcols <- mapmisc::colourScale(
    rast,
    breaks = brk,
    style = "fixed",
    col = "Spectral",
    rev = TRUE,
    dec = -log10(0.05)
  )
  
  if (saveit != "") pdf(saveit,width = 7,height = 7)
  mapmisc::map.new(rast)
  plot(rast,
       col = predcols$col,
       breaks = predcols$breaks,
       legend=FALSE, add=TRUE)
  plot(MGBorder, add=TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
  plot(MGBorderouter,add = TRUE)
  mapmisc::legendBreaks('right', predcols, cex=1.5, bty='n',inset = inset)
  if (saveit != "") dev.off()
}

# Get the mappings from inside their predict functions
tmp_shp <- dis_data$polygon_shapefile
tmp_shp@data <- data.frame(area_id = factor(dis_data$polygon_data$area_id))
shapefile_raster <- raster::rasterize(tmp_shp, dis_data$covariate_rasters, 
                                      field = "area_id")
shapefile_ids <- raster::unique(shapefile_raster)
iid_objects <- list(shapefile_raster = shapefile_raster, 
                    shapefile_ids = shapefile_ids)

iid_ras <- iid_objects$shapefile_raster

coords <- dis_data$coordsForPrediction
Amatrix <- disaggregation:::getAmatrix(dis_data$mesh,coords)

make_eta <- function(W) {
  # W is a column of the matrix of samples
  # It represents a sample from the approximation to the posterior of the parameter
  # vector ordered according to TMB
  # intercept and 3 slopes, 109 polygon-level iid effects, 531 mesh-level random effects
  eta <- 0
  # Covariate effect
  eta <- eta + W[1] + covariate_stack[[1]] * W[2] + covariate_stack[[2]] * W[3] + covariate_stack[[3]] * W[4]
  # eta is now a raster
  # iid effect
  # iideffect <- W[5:113]
  # for (i in seq_along(iideffect)) {
  #   iid_ras@data@values[which(iid_objects$shapefile_raster@data@values == 
  #                               iid_objects$shapefile_ids[i])] <- iideffect[i]
  # }
  # eta <- eta + iid_ras
  # spatial field
  nodemean <- W[114:644]
  field <- Amatrix %*% nodemean
  field_ras <- raster::rasterFromXYZ(cbind(coords,field))
  eta <- eta + field_ras
  eta
}

system.time(ee <- apply(samps$samps,2,make_eta)) # 100 samps: 280 secs
eb <- brick(ee)
predmean <- mean(exp(eb))
predexceedence <- mean(eb > log(.2))

uu <- samps$samps[114:644, ]
uest <- apply(uu,1,mean)
field <- Amatrix %*% uest
uest <- raster::rasterFromXYZ(cbind(coords,field))
covcontrib <- mean(eb) - uest

## Their fit ##
obj$fn(opt$mode) # Eval the random effects at the mode, internally
theirfitW <- obj$env$last.par
theirfiteta <- make_eta(theirfitW)

## Load and process MCMC results ##
mcmcenv <- new.env()
load(mcmcpath,envir = mcmcenv)

mcmcsamples <- as.data.frame(mcmcenv$mcmc_out)

mcmcest <- apply(mcmcsamples,2,mean)

thetaest <- mcmcest[c("iideffect_log_tau","log_sigma","log_rho")]

# Make etas for only 100 samples, for time
set.seed(463279)
samplestodo <- sample(nrow(mcmcsamples),100,FALSE)

West <- mcmcsamples[samplestodo,!(names(mcmcest) %in% c("iideffect_log_tau","log_sigma","log_rho"))]

mcmc_eta <- brick(apply(West,1,make_eta))

mcmc_plotraster <- mean(exp(mcmc_eta))
mcmcexceedence <- mean(mcmc_eta > log(.2))
mcmc_U <- mcmcsamples[samplestodo,117:647]
mcmc_U <- apply(mcmc_U,2,mean)
field <- Amatrix %*% mcmc_U
mcmc_U <- raster::rasterFromXYZ(cbind(coords,field))


## Make the plots ##

plotmap(predmean,seq(0,.45,by=.05),saveit = paste0(figurepath,"predmean-",datestamp,"-",version,".pdf"),inset=.2)
plotmap(predexceedence,seq(0,.95,by=.1),saveit = paste0(figurepath,"predexceedence-",datestamp,"-",version,".pdf"),inset=.2)
plotmap(uest,seq(-2.5,2,by=.5),saveit = paste0(figurepath,"uest-",datestamp,"-",version,".pdf"),inset=.2)

plotmap(mcmc_plotraster,seq(0,.45,by=.05),saveit = paste0(figurepath,"mcmcfit-",datestamp,"-",version,".pdf"),inset=.2)
plotmap(mcmcexceedence,seq(0,.95,by=.1),saveit = paste0(figurepath,"mcmcexceedence-",datestamp,"-",version,".pdf"),inset=.2)
plotmap(mcmc_U,seq(-2.5,2,by=.5),saveit = paste0(figurepath,"mcmcU-",datestamp,"-",version,".pdf"),inset=.2)

# plotmap(population_raster,round(quantile(values(population_raster),probs = (0:10)/10,na.rm=TRUE),-2),saveit = paste0(figurepath,"populationraster-",datestamp,"-",version,".pdf"),inset=.1)


# Polygons and covariates
plotmap(covariate_stack[[1]],seq(-1.2,4.6,by=.5),saveit = paste0(figurepath,"elevation-",datestamp,"-",version,".pdf"),inset=.1)
plotmap(covariate_stack[[2]],seq(-2.8,3.2,by=.5),saveit = paste0(figurepath,"vegetation-",datestamp,"-",version,".pdf"),inset=.1)
plotmap(covariate_stack[[3]],seq(-3.4,2.1,by=.5),saveit = paste0(figurepath,"lst-",datestamp,"-",version,".pdf"),inset=.1)

malariacols <- mapmisc::colourScale(
  dis_data$polygon_data$response,
  breaks = c(200,500,1000,2000,5000,10000,20000,30000,60000),
  dec=-1,
  # style = "equal", transform='log',
  style =  'fixed',
  col = "Oranges"
)
pdf(paste0(figurepath,"polygoncounts-",datestamp,"-",version,".pdf"),width = 7,height = 7)
mapmisc::map.new(dis_data$polygon_shapefile)
plot(dis_data$polygon_shapefile,col = malariacols$plot, add=TRUE, border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(MGBorderouter,add = TRUE)
mapmisc::legendBreaks('topright',malariacols, cex=1.5,bty = 'n',inset=.1)
mapmisc::scaleBar(dis_data$polygon_shapefile, 'topleft', bty='n', cex=2)
dev.off()


# Population
popinpoly <- raster::extract(population_raster,dis_data$polygon_shapefile)
popsuminpoly <- reduce(lapply(popinpoly,sum),c)

options(scipen = 10)
popcols <- mapmisc::colourScale(
  popsuminpoly,
  breaks = round(quantile(reduce(popsuminpoly,c),(0:9)/9,na.rm=TRUE),-4),
  dec=-1,
  # style = "equal", transform='log',
  style =  'fixed',
  col = "Oranges"
)
pdf(paste0(figurepath,"popcounts-",datestamp,"-",version,".pdf"),width = 7,height = 7)
mapmisc::map.new(dis_data$polygon_shapefile)
plot(dis_data$polygon_shapefile,col = popcols$plot, add=TRUE, border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(MGBorderouter,add = TRUE)
mapmisc::legendBreaks('topright',popcols, cex=1.5,bty = 'n',inset=.1)
dev.off()

poprastercols <- mapmisc::colourScale(
  population_raster,
  breaks = round(quantile(values(population_raster),(0:9)/9,na.rm=TRUE),-2),
  dec=-1,
  # style = "equal", transform='log',
  style =  'fixed',
  col = "Oranges"
)

plotraster <- mask(population_raster,dis_data$polygon_shapefile)
pdf(paste0(figurepath,"populationraster-",datestamp,"-",version,".pdf"),width = 7,height = 7)
mapmisc::map.new(plotraster)
plot(plotraster,col = poprastercols$col,
     breaks = poprastercols$breaks, add=TRUE, legend = FALSE, lwd=0.5)
plot(MGBorder, add=TRUE, lwd=0.5)
plot(MGBorderouter,add = TRUE)
mapmisc::legendBreaks('topright',poprastercols, cex=1.5,bty = 'n',inset=.1)
dev.off()

# Smoothing param posteriors
smoothingparamsamples <- as_tibble(mcmcsamples[ ,c("iideffect_log_tau","log_sigma","log_rho")])

smoothingparam_tmbest <- sdsummary

tauplotmcmc <- smoothingparamsamples %>%
  ggplot(aes(x = exp(-.5*iideffect_log_tau))) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 100,colour = 'transparent',fill = 'lightgray') +
  # ELGM approx
  geom_line(data = taupdf,mapping = aes(x = transparam,y = pdf_transparam)) +
  # Gaussian approx
  stat_function(fun = dlnorm,
                args = list(
                  meanlog = -.5*smoothingparam_tmbest['iideffect_log_tau','Estimate'],
                  sd = .5*smoothingparam_tmbest['iideffect_log_tau','Std. Error']
                ),
                linetype = 'dashed'
  ) +
  stat_function(fun = sigmaprior,args = list(u = tau_u,alpha = tau_alpha),linetype = 'dotdash') +
  coord_cartesian(xlim = c(0.25,1)) +
  labs(x = expression(1/sqrt(tau)),y = '') +
  theme(text = element_text(size = NLPLOTTEXTSIZE))

sigmaplotmcmc <- smoothingparamsamples %>%
  ggplot(aes(x = exp(log_sigma))) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 100,colour = 'transparent',fill = 'lightgray') +
  # ELGM approx
  geom_line(data = sigmapdf,mapping = aes(x = transparam,y = pdf_transparam)) +
  # Gaussian approx
  stat_function(fun = dlnorm,
                args = list(
                  meanlog = smoothingparam_tmbest['log_sigma','Estimate'],
                  sd = smoothingparam_tmbest['log_sigma','Std. Error']
                ),
                linetype = 'dashed'
  ) +
  stat_function(fun = sigmaprior,args = list(u = sigma_u,alpha = sigma_alpha),linetype = 'dotdash') +
  coord_cartesian(xlim = c(0.25,2)) +
  labs(x = expression(sigma),y = '') +
  theme(text = element_text(size = NLPLOTTEXTSIZE))

rhoplotmcmc <- smoothingparamsamples %>%
  filter(exp(log_rho) < 6) %>%
  ggplot(aes(x = exp(log_rho))) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),bins = 100,colour = 'transparent',fill = 'lightgray') +
  # ELGM approx
  geom_line(data = filter(rhopdf,transparam < 6),mapping = aes(x = transparam,y = pdf_transparam)) +
  # Gaussian approx
  stat_function(fun = dlnorm,
                args = list(
                  meanlog = smoothingparam_tmbest['log_rho','Estimate'],
                  sd = smoothingparam_tmbest['log_rho','Std. Error']
                ),
                linetype = 'dashed'
  ) +
  stat_function(fun = rhoprior,args = list(u = rho_u,alpha = rho_alpha),linetype = 'dotdash') +
  labs(x = expression(rho),y = '') +
  theme(text = element_text(size = NLPLOTTEXTSIZE))

WIDTH <- 4
HEIGHT <- 2.5

ggsave(
  paste0(figurepath,"tauplot-",datestamp,"-",version,".pdf"),
  tauplotmcmc,
  width = WIDTH,height = HEIGHT
)

ggsave(
  paste0(figurepath,"sigmaplot-",datestamp,"-",version,".pdf"),
  sigmaplotmcmc,
  width = WIDTH,height = HEIGHT
)

ggsave(
  paste0(figurepath,"rhoplot-",datestamp,"-",version,".pdf"),
  rhoplotmcmc,
  width = WIDTH,height = HEIGHT
)

# tauplotmcmc / sigmaplotmcmc / rhoplotmcmc

## Summary of regression coefficients ##

beta_idx <- c(1:4)
nonlinear_idx <- (114:116)

mcmcbetamean <- apply(mcmcsamples[ ,beta_idx],2,mean)
mcmcbetasd <- apply(mcmcsamples[ ,beta_idx],2,sd)
mcmcbeta2.5 <- apply(mcmcsamples[ ,beta_idx],2,quantile,probs = .025)
mcmcbeta97.5 <- apply(mcmcsamples[ ,beta_idx],2,quantile,probs = .975)

mcmcnonlinearsamps <- cbind(
  exp(-.5*mcmcsamples[ ,nonlinear_idx[1]]),
  exp(mcmcsamples[ ,nonlinear_idx[2]]),
  exp(mcmcsamples[ ,nonlinear_idx[3]])
)

mcmcnonlinearmean <- apply(mcmcnonlinearsamps,2,mean)
mcmcnonlinearsd <- apply(mcmcnonlinearsamps,2,sd)
mcmcnonlinear2.5 <- apply(mcmcnonlinearsamps,2,quantile,probs = .025)
mcmcnonlinear97.5 <- apply(mcmcnonlinearsamps,2,quantile,probs = .975)


# Use the same number of iid post samps
set.seed(438092)
samps_forbeta <- sample_marginal(aghqmodel,nrow(mcmcsamples))
betasamps <- samps_forbeta$samps[1:4, ]

elgmbetamean <- apply(betasamps,1,mean)
elgmbetasd <- apply(betasamps,1,sd)
elgmbeta2.5 <- apply(betasamps,1,quantile,probs = .025)
elgmbeta97.5 <- apply(betasamps,1,quantile,probs = .975)

elgmnonlinearmean <- compute_moment(aghqmodel$normalized_posterior,function(x) c(exp(-.5*x[1]),exp(x[2:3])))
elgmnonlinearsd <- sqrt(compute_moment(
  aghqmodel$normalized_posterior,
  function(x) c( (exp(-.5*x[1]) - elgmnonlinearmean[1])^2,
                 (exp(x[2:3]) - elgmnonlinearmean[2:3])^2
  )
))
elgmquantstau <- sort(exp(-.5*compute_quantiles(aghqmodel$marginals[[1]])))
elgmquantssigma <- exp(compute_quantiles(aghqmodel$marginals[[2]]))
elgmquantsrho <- exp(compute_quantiles(aghqmodel$marginals[[3]]))

elgmquants2.5 <- c(elgmquantstau[1],elgmquantssigma[1],elgmquantsrho[1])
elgmquants97.5 <- c(elgmquantstau[2],elgmquantssigma[2],elgmquantsrho[2])


coeftable <- tibble(
  meanelgm = c(elgmbetamean,elgmnonlinearmean),
  meanmcmc = c(mcmcbetamean,mcmcnonlinearmean),
  sdelgm = c(elgmbetasd,elgmnonlinearsd),
  sdmcmc = c(mcmcbetasd,mcmcnonlinearsd),
  q2.5elgm = c(elgmbeta2.5,elgmquants2.5),
  q2.5mcmc = c(mcmcbeta2.5,mcmcnonlinear2.5),
  q97.5elgm = c(elgmbeta97.5,elgmquants97.5),
  q97.5cmc = c(mcmcbeta97.5,mcmcnonlinear97.5)
)

readr::write_csv(coeftable,paste0(figurepath,"coefficienttable-",datestamp,"-",version,".csv"))

# Help out with the latex...
knitr::kable(
  coeftable,
  format = "latex",
  digits = 3
)



