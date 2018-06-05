## This model calibrates the bird ecology model from Thiele et al (2014) using MCMC on a model for the population growth rate

# import libraries
library(reticulate) # provides interface with Python

## ======================================================
## compute log likelihood for AR1 errors centered at 0
## X(t) = rho*X(t-1) + eps(t), eps(t) ~ N(0, sigma)
## r are the residuals

logl.ar1 <- function(r, sigma, rho) {
  
  n <- length(r)
  logl <- 0
  if (n > 1) {
    w <- r[2:n] - rho*r[1:(n-1)] # whiten the residuals
    logl <- logl + sum(dnorm(w, sd=sigma, log=TRUE)
  }
  
  logl
}

## =====================================================
## run model with parameters pars
## pass model module as an R object (imported using reticulate)

run_model <- function(scp, sup, ssp, b.mod, n_years, n_cells) {
  # extract parameters

  # construct model instance
  m <- b.mod$BirdModel(scout_prob=scp, surv_prob=sup, scout_surv_prob=ssp, num_years=as.integer(n_years), num_territory=as.integer(n_cells), seed=as.integer(0), model_queries=list(pop=b.mod$get_pop))
  m$run_model
  
  # return population series
  m$query_out$model_vals$pop
}

## =====================================================
## compute log likelihood for the overall model
## pass model module as an R object (using reticulate),
##    model parameters, model parameter names, and aux should be
##    a named list consisting of the number of years and number of cells

log.lik <- function(pars, parnames, obs, b.mod, aux) {
  
  # extract parameters
  # regression parameters
  a <- pars[match('a', parnames)]
  b <- pars[match('b', parnames)]
  # model parameters
  scp <- pars[match('scout_prob', parnames)]
  sup <- pars[match('surv_prob', parnames)]
  ssp <- pars[match('scout_surv_prob', parnames)]
  # AR1 parameters
  rho <- pars[match('rho', parnames)]
  sigma <- pars[match('sigma', parnames)]
  
  # auxiliary model parameters
  n_yr <- aux[['n_years']
  n_terr <- aux[['n_cells']
  
  # run model with sampled model parameters
  pop <- run_model(scp, sup, ssp, b.mod, n_yr, n_terr)
  # compute growth rate for model output
  x <- (pop[2:n_yr] - pop[1:(n_yr-1)])/pop[1:(n_yr-1)]
  r <- obs - (a + b*x) # compute residuals
  
  # return log likelihood
  logl.ar1(r, sigma, rho)
}

## =========================================================
## compute prior log-densities for sampled values
## prior families and hyperparameters should be passed as a
##    named list

log.pri <- function(pars, parnames, priors) {
  
  lpri <- 0
  # loop over parameters
  for (i in 1:length(parnames) {
    name <- parnames[i] # get parameter name
    # if parameters are outside of the range with positive support,
    #   set the log-prior to -Inf, otherwise compute prior value
    if (priors[[name]][['type']] == 'beta') {
      if ((pars[i] <= 0) || (pars[i] >= 1)) {
        lpri <- -Inf
      } else {
        lpri <- lpri + dbeta(pars[i], shape1=priors[[name]][['a']], shape2=priors[[name]][['b']], log=TRUE)
      }
    } else if (priors[[name]][['type']] == 'half.cauchy') {
      if (pars[i] <= 0) {
        lpri <- -Inf
      } else {
        p <- dcauchy(pars[i], location=priors[[name]][['loc']], scale=priors[[name]][['scale']], log=TRUE)
        lpri <- lpri + 2*p # we use half-Cauchy instead of Cauchy
      }
    } else if (priors[[name]][['type']] == 'normal') {
      lpri <- lpri + dnorm(pars[i], mean=priors[[name]][['mean']], sd=priors[[name]][['sd']], log=TRUE)
    } else if (priors[[name]][['type']] == 'uniform') {
      if ((pars[i] <= priors[[name]][['bound.low']]) || (pars[i] >= priors[[name]][['bound.high']])) {
        lpri <- -Inf
      } else {
        lpri <- lpri + dunif(pars[i], min=priors[[name]][['bound.low']], max=priors[[name]][['bound.high']], log=TRUE)
      }
    }
    
  }
  # return log-prior
  lpri
}

## =====================================================
## compute posterior log-density for sampled values

log.post <- function(pars, parnames, priors, obs, b.mod, aux) {
  
  # compute log-prior
  lpri <- log.pri(pars, parnames, priors)
  
  # if log-prior is -Inf, no reason to run the model; otherwise,
  # compute the log-likelihood
  if (is.finite(lpri)) {
    llik <- log.lik(pars, parnames, obs, b.mod, aux)
    lpost <- llik + lpri
  } else {
    lpost <- -Inf
  }
  # return log-posterior
  lpost
}