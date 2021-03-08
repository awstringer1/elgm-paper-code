### Opioid treatment data, town-level random intercept ###
# This script implements data cleaning and defines functions for
# model fitting, supporting the Opioid treatment data example


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
  
  # glimpse(opioid_clean)
  
  # Save the data to rds for faster loading
  readr::write_rds(opioid_clean,savepath,compress = "none")
  
  # Note that removing the missing CBSA values dropped 2,545,961 records, or 26%
}
## Function to fit model ----
# Write a function that fits the model using our approach with or without
# epslion
# This will be run on subsets of the data of various sizes and compared to INLA

fit_model <- function(dat,
                      withepsilon = FALSE,
                      k = 3) {
  # dat: dataset to use. This will be the opioid_clean data created above, 
  # but subsetted according to the size of comparison being run
  # withepsilon: logical, include the epsilon or not? Not used currently.
  # k: number of quadrature points to use
  
  ## Prep design matrix ##
  us <- unique(dat$state)
  ut <- unique(dat$town)
  lss <- length(dat$state)
  ltt <- length(dat$town)

  Zmat <- cbind(
    # Random effects
    cbind(
      Diagonal(n = lss)[match(dat$state,us),order(us)],
      Diagonal(n = ltt)[match(dat$town,ut),order(ut)]
    ),
    # Fixed effects
    sparse.model.matrix(completed ~ gender + race + livingarrangement,data = dat)
  )
  
  # Dimensions
  Udim <- c(length(us),length(ut)) # Dimension of random effects
  Wd <- ncol(Zmat) # Total dimension
  p <- Wd - sum(Udim) # Number of fixed effects
  n <- nrow(dat)
  if (withepsilon) Wd <- Wd + n # Add the epsilon
  
  ## Priors ----
  # Prior on theta = -2log(sigma), sigma = exp(-theta/2)
  log_prior_theta <- function(theta,prior_alpha = c(.5,.5),prior_u = c(.5,.5)) {
    # Implement the PC prior directly.
    # P(sigma > u) = alpha.
    # See inla.doc("pc.prec")
    lambda <- -log(prior_alpha)/prior_u
    sum(log(lambda/2) - lambda * exp(-theta/2) - theta/2)
  }
  
  # Prior on beta is gaussian with precision...
  prior_prec_beta <- .001

  # Q-matrix, prior precision matrix of Gaussians
  if (!withepsilon) {
    Q_matrix <- function(theta) {
      Diagonal(n = Wd,
               x = c(
                 rep(exp(theta[1]),Udim[1]),
                 rep(exp(theta[2]),Udim[2]),
                 rep(prior_prec_beta,p)
               )
      )
    }
  } else {
    tau <- exp(12) # Default in the R-INLA software
    Q_matrix <- function(theta) {
      QQ <- Diagonal(n = sum(Udim) + p,
                     x = c(
                       rep(exp(theta[1]),Udim[1]),
                       rep(exp(theta[2]),Udim[2]),
                       rep(prior_prec_beta,p)
                     )
               )
      tau * rbind(
        cbind(Diagonal(n = n),-Zmat),
        cbind(-t(Zmat),(1/tau) * QQ + crossprod(Zmat))
      )
    }
  }

  
  log_prior_W <- function(W,theta) {
    Q <- Q_matrix(theta)
    dt <- as.numeric(determinant(Q,logarithm = TRUE)$modulus) # Diagonal matrix, so very fast
    (1/2) * dt - (1/2) * as.numeric(crossprod(W,crossprod(Q,W)))
  }
  
  ## Likelihood ----
  if (!withepsilon) {
    make_eta <- function(W) as.numeric(Zmat %*% W)
  } else {
    make_eta <- function(W) W[1:n]
  }
  
  log_likelihood_binomial <- function(W) {
    eta <- make_eta(W)
    sum(dbinom(dat$completed,1,prob = (exp(eta)/(1 + exp(eta))),log = TRUE))
  }
  
  grad_loglik_binomial <- function(W) {
    eta <- make_eta(W)
    dat$completed - exp(eta)/(1 + exp(eta))
  }
  
  # Negative hessian of loglik wrt eta
  C_matrix <- function(W) {
    eta <- make_eta(W)
    p <- exp(eta)/(1+exp(eta))
    Diagonal(length(eta),p * (1 - p))
  }
  
  ## Posterior ----

  log_posterior_W <- function(W,theta) {
    log_prior_W(W,theta) + log_likelihood_binomial(W)
  }
  
  if (!withepsilon) {
    grad_log_posterior_W <- function(W,theta) {
      as.numeric(-Q_matrix(theta) %*% cbind(W) + crossprod(Zmat,grad_loglik_binomial(W)))
    }
    
    # POSITIVE hessian of log-post wrt W
    H_matrix <- function(W,theta) {
      Ceta <- C_matrix(W)
      CW <- crossprod(Zmat,crossprod(Ceta,Zmat))
      
      -1 * as(Q_matrix(theta) + CW,"dgCMatrix")
    }
  } else {
    # If epsilon, now gradient and hessian are wrt eta and then zero-padded
    grad_log_posterior_W <- function(W,theta) as.numeric(-Q_matrix(theta) %*% cbind(W) + c(grad_loglik_binomial(W),rep(0,Wd - n)))
    H_matrix <- function(W,theta) {
      Ceta <- C_matrix(W)
      -1 * (Q_matrix(theta) + bdiag(Ceta,Diagonal(n = Wd - n,x = 0)))
    }
  }
  
  ## AGHQ ----
  startingvalues <- list(W = rep(0,Wd),theta = c(0,0))
  ff <- list(
    fn = log_posterior_W,
    gr = grad_log_posterior_W,
    he = H_matrix
  )
  
  # Fit the model
  modelobject <- aghq::marginal_laplace(ff,k,startingvalues)
  # Include time taken for posterior summaries and samples in the stated runtime
  summary(modelobject)
  aghq::sample_marginal(modelobject,1e03) # 1000 posterior samples
  
  return(modelobject) # Prevent returning the big set of samples
}

## INLA ----

fit_inla_multi <- function(dat,threads) {
  n <- nrow(dat)
  inlamod_multi <- tryCatch(inla(
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
    num.threads = threads,
    blas.num.threads = threads
  ),
  error = function(e) {print(paste0("Error when fitting INLA multi threaded, n = ",n,"\n")); return(e)}
  )
  inlamod_multi
}

fit_inla_single <- function(dat) {
  n <- nrow(dat)
  inlamod_single <- tryCatch(inla(
    completed ~ gender + race + livingarrangement + 
      f(state,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))) +
      f(town,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))),
    data = dat,
    family = 'binomial',
    control.compute = list(
      openmp.strategy = 'huge',
      smtp = 'pardiso'
    ),
    control.inla = list(
      strategy = 'gaussian',
      int.strategy = 'ccd'
    ),
    num.threads = 1,
    blas.num.threads = 1
  ),
  error = function(e) {print(paste0("Error when fitting INLA single threaded, n = ",n,"\n")); return(e)}
  )
  inlamod_single
}
