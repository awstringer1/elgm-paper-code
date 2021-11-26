### Astro example ###
# Reproduce example 6.3
# Data obtained from aghq package
# Results saved to temp storage
# This script should be sourcable without modification


## BEGIN SETUP ##

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'Matrix',
  'aghq',
  'TMB'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

## Setup ----

globalpath <- tempdir()
figurepath <- file.path(globalpath,"figures")
if (!dir.exists(figurepath)) dir.create(figurepath)
datestamp <- "20210211"
version <- "v1"

# Precompile TMB
precompile()
# Compile and load the TMB template. Source code included in the aghq package
tmbpath <- tempdir()
stopifnot(file.copy(system.file('extsrc/04_astro.cpp',package='aghq'),tmbpath))
compile(file.path(tmbpath,"04_astro.cpp"))
dyn.load(dynlib(file.path(tmbpath,"04_astro")))

## Load data ----
data("gcdatalist",package = 'aghq')

## Set up model ----

parambounds <- list(
  Psi0 = c(1,200),
  gamma = c(.3,.7),
  alpha = c(3.0,3.7),
  beta = c(-.5,1)
)

get_psi0 <- function(theta) {
  (parambounds$Psi0[2] - parambounds$Psi0[1]) * exp(-exp(theta)) + parambounds$Psi0[1]
}
get_theta1 <- function(Psi0) log(-log( (Psi0 - parambounds$Psi0[1]) / (parambounds$Psi0[2] - parambounds$Psi0[1]) ))

get_gamma <- function(theta) {
  (parambounds$gamma[2] - parambounds$gamma[1]) * exp(-exp(theta)) + parambounds$gamma[1]
}
get_theta2 <- function(gamma) log(-log( (gamma - parambounds$gamma[1]) / (parambounds$gamma[2] - parambounds$gamma[1]) ))

get_alpha <- function(theta) {
  exp(theta) + parambounds$alpha[1]
}
get_theta3 <- function(alpha) log(alpha - parambounds$alpha[1])

get_beta <- function(theta) {
  (parambounds$beta[2] - parambounds$beta[1]) * exp(-exp(theta)) + parambounds$beta[1]
}
get_theta4 <- function(beta) log(-log( (beta - parambounds$beta[1]) / (parambounds$beta[2] - parambounds$beta[1]) ))

get_theta <- function(params) {
  c(
    get_theta1(params[1]),
    get_theta2(params[2]),
    get_theta3(params[3]),
    get_theta4(params[4])
  )
}

get_params <- function(theta) {
  c(
    get_psi0(theta[1]),
    get_gamma(theta[2]),
    get_alpha(theta[3]),
    get_beta(theta[4])
  )
}


# Function to create the AD object for fixed theta
get_ad_object <- function(theta,constraintfun = FALSE) {
  # Organize the data and call MakeADFun for the given theta
  # Confusing terminology; in this context, "theta" is "data", at least
  # from the perspective of TMB
  # 
  # If constraintfun = TRUE then the ADReport = TRUE flag is used and
  # an AD object corresponding to the nonlinear constraints is returned instead
  w <- matrix(0,nrow = nrow(gcdatalist$y),ncol = 2) # Don't include measurement errors for PM right now
  
  paramstart <- get_params(theta)
  datlist <- gcdatalist
  
  datlist$p <- paramstart[1]
  datlist$g <- paramstart[2]
  datlist$a <- paramstart[3]
  datlist$b <- paramstart[4]
  
  ff <- MakeADFun(data = datlist,parameters = list(W = w),DLL = "04_astro",ADreport = constraintfun,silent = TRUE)
  
  ff
}

make_w_matrix <- function(W) {
  n <- length(W)/2
  
  idx1 <- 1:n
  idx2 <- idx1 + n

  cbind(W[idx1],W[idx2]) # Don't include measurement errors for PM right now
}

log_likelihood <- function(W,ff,theta) {
  # W: vector of length 4n containing measurement errors 1:n for
  # Rgc, Vlos, PMra and PMdec, in that order
  # ff: an AD object returned by get_ad_obj
  # theta: nonlinear parameters
  # 
  # Add on the normconst (does not depend on W, but used in theta part)
  # For some reason I couldn't get that to work in TMB
  
  params <- get_params(theta)
  p <- params[1]
  g <- params[2]
  a <- params[3]
  b <- params[4]
  
  # n <- length(W) / 4
  n <- length(W) / 2
  
  ww <- make_w_matrix(W)
  ll <- ff$fn(ww)
  normconst <- (2*b/g - a/g)*log(p) - log( pi*sqrt(pi)*2^( -b+3/2 ) ) - lgamma(x = (1-b)) +
    lgamma( x = (a/g - 2*b/g + 1) ) - lgamma( x = ( b*(g-2)/g + a/g - 1/2 ) )
  ll + n * normconst
  # ll
}

grad_log_likelihood <- function(W,ff) {
  # W: vector of length 4n containing measurement errors 1:n for
  # Rgc, Vlos, PMra and PMdec, in that order
  # ff: an AD object returned by get_ad_obj
  
  ww <- make_w_matrix(W)
  as.numeric(ff$gr(ww))
}

hess_log_likelihood <- function(W,ff) {
  # W: vector of length 4n containing measurement errors 1:n for
  # Rgc, Vlos, PMra and PMdec, in that order
  # ff: an AD object returned by get_ad_obj
  
  ww <- make_w_matrix(W)
  as(ff$he(ww),"sparseMatrix")
}

## Prior + Posterior ###

# The Q-matrix in this model is FIXED, i.e. it does not depend on the hyperparameters
# Keep it as a "function" to avoid confusion with previous code. But here, theta does nothing.
Q_matrix <- function(theta) {
  sd = c(
    gcdata$eRgc,
    gcdata$eVlos # Don't include measurement errors for PM right now
    # gcdata$ePMra,
    # gcdata$ePMdec
  )
  prec = (1/sd)^2
  # Diagonal(n = 4 * nrow(gcdata),x = prec)
  Diagonal(n = 2 * nrow(gcdata),x = prec)
}
# For places in the code which use precomputed "Q"... precompute Q
Q <- Q_matrix(0)
Qdet <- as.numeric(determinant(Q,logarithm = TRUE)$modulus)
# but it takes less than a millisecond so it doesn't matter too much

log_prior_W <- function(W,theta) {
  -(1/2) * as.numeric(crossprod(W,crossprod(Q,W))) + (1/2)*Qdet
}

# Prior for theta-- long construction, compared to other models
# Most of it is done above
# All priors
theta1logprior <- function(theta) dexp(exp(theta),1,log = TRUE) + theta
theta2logprior <- function(theta) dexp(exp(theta),1,log = TRUE) + theta
theta4logprior <- function(theta) dexp(exp(theta),1,log = TRUE) + theta
# log-gamma for theta3 = log(alpha - 3)
theta3logprior <- function(theta) theta + dgamma(exp(theta),shape = 1,rate = 4.6,log = TRUE)
log_prior_theta <- function(theta) {
  # theta = c(theta1,theta2,theta3,theta4)
  theta1logprior(theta[1]) +
    theta2logprior(theta[2]) +
    theta3logprior(theta[3]) +
    theta4logprior(theta[4])
}

log_posterior_W <- function(W,theta) {
  ff <- get_ad_object(theta)
  log_prior_W(W,theta) + log_likelihood(W,ff,theta) + log_prior_theta(theta)
}

grad_log_posterior_W <- function(W,theta) {
  ff <- get_ad_object(theta)
  as.numeric(-Q %*% W) + grad_log_likelihood(W,ff)
}

# POSITIVE hessian
H_matrix <- function(W,theta) {
  ff <- get_ad_object(theta)
  -Q + hess_log_likelihood(W,ff)
}

## Fit the model ----

ff <- list(
  fn = log_posterior_W,
  gr = grad_log_posterior_W,
  he = H_matrix
)
cat("Fitting model.\n")
aghqmodel <- marginal_laplace(
  ff,
  k = 3,
  startingvalue = list(W = rep(0,2*nrow(gcdatalist$y)),theta = get_theta(c(26.000,0.375,3.040,0.400)))
)

## Posterior summaries ----

cat("Computing summaries.\n")

PLOTTEXTSIZE <- 28

Psi0prior <- function(Psi0) dunif(Psi0,parambounds$Psi0[1],parambounds$Psi0[2],log = FALSE)
gammaprior <- function(gamma) dunif(gamma,parambounds$gamma[1],parambounds$gamma[2],log = FALSE)
alphaprior <- function(alpha) dgamma(alpha - parambounds$alpha[1],shape = 1,rate = 4.6,log = FALSE)
betaprior <- function(beta) dunif(beta,parambounds$beta[1],parambounds$beta[2],log = FALSE)


psi0post <- compute_pdf_and_cdf(
  aghqmodel$marginals[[1]],
  list(totheta = get_theta1,fromtheta = get_psi0)
)

# Add a little buffer for numerical stability
get_theta2 <- function(gamma) 
  log(-log( (gamma - parambounds$gamma[1] + 1e-03) / 
              (parambounds$gamma[2] - parambounds$gamma[1] + 1e-03) ))

gammapost <- compute_pdf_and_cdf(
  aghqmodel$marginals[[2]],
  list(totheta = get_theta2,fromtheta = get_gamma)
)

get_theta3 <- function(alpha) 
  log(alpha - parambounds$alpha[1] + 1e-03)

alphapost <- compute_pdf_and_cdf(
  aghqmodel$marginals[[3]],
  list(totheta = get_theta3,fromtheta = get_alpha)
)

betapost <- compute_pdf_and_cdf(
  aghqmodel$marginals[[4]],
  list(totheta = get_theta4,fromtheta = get_beta)
)


psi0_postplot <- psi0post %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_line(size = 1) +
  stat_function(fun = Psi0prior,linetype = "dashed") +
  labs(title = "",x = expression(Psi),y = "") +
  # labs(title = "",x = "",y = "") +
  scale_x_continuous(breaks = seq(24,42,by = 2)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

gamma_postplot <- gammapost %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_line(size = 1) +
  stat_function(fun = gammaprior,linetype = "dashed") +
  labs(title = "",x = expression(gamma),y = "") +
  # labs(title = "",x = "",y = "") +
  scale_x_continuous(breaks = seq(.3,.45,by = .02)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

alpha_postplot <- alphapost %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_line(size = 1) +
  stat_function(fun = alphaprior,linetype = "dashed") +
  labs(title = "",x = expression(alpha),y = "") +
  # labs(title = "",x = "",y = "") +
  scale_x_continuous(breaks = seq(3,3.3,by = .02)) +
  theme(text = element_text(size = PLOTTEXTSIZE)) +
  coord_cartesian(xlim = c(3,3.1))

beta_postplot <- betapost %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_line(size = 1) +
  stat_function(fun = betaprior,linetype = "dashed") +
  labs(title = "",x = expression(beta),y = "") +
  # labs(title = "",x = "",y = "") +
  scale_x_continuous(breaks = seq(-.3,.4,by = .2)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

## Cumulative Mass ----

Mr <- function(r,theta) {
  p = get_psi0(theta[1])
  g = get_gamma(theta[2])
  
  # Manual unit conversion into "mass of one trillion suns" (so awesome)
  g*p*r^(1-g) * 2.325e09 * 1e-12
}

rtodo <- 1:150
Mrout <- numeric(length(rtodo))
Mrsdout <- numeric(length(rtodo))
for (rr in 1:length(rtodo)) {
  r <- rtodo[rr]
  
  Mrout[rr] <- compute_moment(
    aghqmodel$normalized_posterior,
    function(x) Mr(r,x)
  )
  Mrsdout[rr] <- sqrt(compute_moment(
    aghqmodel$normalized_posterior,
    function(x) (Mr(r,x) - Mrout[rr])^2
  ))
}

massplot <- tibble(
  r = rtodo,
  Mr = Mrout,
  Mrsd = Mrsdout
) %>%
  ggplot(aes(x = r)) +
  theme_classic() +
  geom_line(aes(y = Mr),size = 1) +
  geom_ribbon(aes(ymin = Mrout - Mrsdout,
                  ymax = Mrout + Mrsdout),
              alpha = .5) +
  geom_ribbon(aes(ymin = Mrout - 2*Mrsdout,
                  ymax = Mrout + 2*Mrsdout),
              alpha = .2) +
  labs(title = "",
       x = "r (kpc)",
       y = bquote('M(r) ('~10^12~M[sun]~')')) +
  scale_x_continuous(breaks = seq(0,150,by = 25)) +
  scale_y_continuous(breaks = seq(0,1,by=.1)) +
  theme(text = element_text(size = PLOTTEXTSIZE))

# Save the plots
WIDTH <- HEIGHT <- 7
ggsave(
  file.path(figurepath,paste0("psi0postplot-",datestamp,"-",version,".pdf")),
  psi0_postplot,
  width = WIDTH,height = HEIGHT
)
ggsave(
  file.path(figurepath,paste0("gammapostplot-",datestamp,"-",version,".pdf")),
  gamma_postplot,
  width = WIDTH,height = HEIGHT
)
ggsave(
  file.path(figurepath,paste0("alphapostplot-",datestamp,"-",version,".pdf")),
  alpha_postplot,
  width = WIDTH,height = HEIGHT
)
ggsave(
  file.path(figurepath,paste0("betapostplot-",datestamp,"-",version,".pdf")),
  beta_postplot,
  width = WIDTH,height = HEIGHT
)
ggsave(
  file.path(figurepath,paste0("masspostplot-",datestamp,"-",version,".pdf")),
  massplot,
  width = WIDTH,height = HEIGHT
)


cat(paste0("Done. You can go to ",figurepath," to see the output. Thanks!\n"))

