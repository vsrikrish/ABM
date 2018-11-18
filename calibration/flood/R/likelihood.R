# contains functions for computing the posterior density of the flood model
library(plyr)
source('utils.R')
##################################################
# ev_mod: computes expected value for likelihood #
#                                                #
# Inputs:                                        #
#   1) model: model func          #
#   2) parnames: vector of parameter names       #
#   3) priors: list with parnames as names       #
#              and each entry is a list with     #
#              values FUN (density function)     #
#              and params (list of parameters)   #
#                                                #
# Output:                                        #
#   1) log-prior density value                   #
##################################################
ev_mod <- function(pars, parnames, type, ffreq) {
  # assign parameters to variables by name
  int <- pars[match('int', parnames)]
  scoef <- pars[match('self_coef', parnames)]
  ncoef <- pars[match('nghd_coef', parnames)]
  fill <- pars[match('fill', parnames)]

  # construct parameter list for model run
  params <- list('int'=int, 'self_coef'=scoef, 'fill_prob'=fill)
  if (type == 'complex') {
    params['nghd_coef'] <- ncoef
  }
  # run model
#  model_out <- run_model(model)(params)
  probs <- vector('list', dim(ffreq)[3])
  if (type == 'simple') {
    for (t in 1:dim(ffreq)[3]) {
      if (t == 1) {
        probs[[t]] <- apply(inv_logit(int+scoef*ffreq[,,1]/10), c(1,2), function(p) {(c(0.99, 0.01) %*% matrix(c(1-p, p, fill, 1-fill), nrow=2,byrow=T))[,1]})
      } else {
        probs[[t]] <- matrix(laply(Map(function(x, y) x %*% y, lapply(alply(probs[[t-1]], c(1,2)), function(p) c(p, 1-p)), alply(inv_logit(int+scoef*ffreq[,,t]/10), c(1, 2), function(p) {matrix(c(1-p, p, fill, 1-fill), nrow=2,byrow=T)})), function(l) l[1]), nrow=dim(ffreq)[1], ncol=dim(ffreq)[2])
      }
    }

  } else {
    # allocate space for transition probabilities and neighborhood storage
    coords <- expand.grid(1:dim(ffreq)[1], 1:dim(ffreq)[2])
    nghd <- vector('list', dim(ffreq)[1]*dim(ffreq)[2])
    for (i in 1:prod(dim(ffreq)[1:2])) {
        ngh <- intersect(which(abs(coords[i,1] - coords[,1]) <= 1),
                          which(abs(coords[i,2] - coords[,2] )<= 1))
        nghd[[i]] <- ngh[!ngh == i]
    }

    for (t in 1:dim(ffreq)[3]) {
      if (t == 1) {
        probs[[t]] <- apply(inv_logit(int+scoef*ffreq[,,1]/10+ncoef*0.01), c(1,2), function(p) {(c(0.99, 0.01) %*% matrix(c(1-p, p, fill, 1-fill), nrow=2,byrow=T))[,1]})
      } else {
        nghd_vac <- 1-matrix(vapply(1:prod(dim(ffreq)[1:2]), function(i) {mean(as.numeric(probs[[t-1]])[nghd[[i]]])}, numeric(1)), nrow=dim(ffreq)[1], ncol=dim(ffreq)[2])
        probs[[t]] <- matrix(laply(Map(function(x, y) x %*% y, lapply(alply(probs[[t-1]], c(1,2)), function(p) c(p, 1-p)), alply(inv_logit(int+scoef*ffreq[,,t]/10+ncoef*nghd_vac), c(1, 2), function(p) {matrix(c(1-p, p, fill, 1-fill), nrow=2,byrow=T)})), function(l) l[1]), nrow=dim(ffreq)[1], ncol=dim(ffreq)[2])
      }
    }

  }

  aperm(laply(probs, identity), c(2, 3, 1))

}

##################################################
# log_prior: computes prior log-value            #
#                                                #
# Inputs:                                        #
#   1) pars: vector of parameter values          #
#   2) parnames: vector of parameter names       #
#   3) priors: list with parnames as names       #
#              and each entry is a list with     #
#              values FUN (density function)     #
#              and params (list of parameters)   #
#                                                #
# Output:                                        #
#   1) log-prior density value                   #
##################################################
log_prior <- function(pars, parnames, priors) {
  lp <- 0 # initialize variable for log prior values
  # loop over parameters
  for (name in parnames) {
    # find parameter value and call prior density function
    par <- pars[match(name, parnames)]
    lp <- lp + do.call(match.fun(priors[[name]]$FUN),
                  c(list(x=par), priors[[name]]$params, list(log=TRUE))
               )
  }
  lp # return log prior value
}


##################################################
# log_lik: compute Poisson log-likelihood based  #
#   on observed quantities                       #
#                                                #
# Inputs:                                        #
#   1) pars: vector of parameter values          #
#   2) parnames: vector of parameter names       #
#   3) model: model function object              #
#   4) dat: vector of observation counts         #
#   5) args: named list of required positional   #
#            arguments, which must include:      #
#       a) n_col: number of columns in grid      #
#       b) n_row: number of rows in grid         #
#       c) mem_length: length of agent memory    #
#       d) seed: random seed for python/numpy    #
#   6) ...: keyword arguments to pass to Python  #
#           API in keyword-value pairs           #
#                                                #
# Output:                                        #
#   1) log-likelihood value                      #
##################################################
log_lik_ind <- function(pars, parnames, ev_fun, dat, type, ffreq)
{
  force(ev_fun)

  probs <- ev_fun(pars, parnames, type, ffreq)

  dbinom(dat, size=1, prob=probs, log=TRUE)  # calculate log-likelihood and return

}

log_lik_agg <- function(pars, parnames, ev_fun, dat, type, ffreq)
{
  force(ev_fun)

  probs <- ev_fun(pars, parnames, type, ffreq)
  exp_vac <- prod(dim(ffreq)[1:2]) - apply(probs, 3, sum)

  dpois(dat, lambda=exp_vac, log=TRUE)
}

##################################################
# neg_log_lik: compute the negative              #
#             log-likelihood (for use with       #
#             optim functions)                   #
#                                                #
# Inputs:                                        #
#   1) pars: vector of parameter values          #
#   2) parnames: vector of parameter names       #
#   3) model: model function object              #
#   4) dat: vector of observation counts         #
#   5) args: named list of required positional   #
#            arguments, which must include:      #
#       a) n_col: number of columns in grid      #
#       b) n_row: number of rows in grid         #
#       c) mem_length: length of agent memory    #
#       d) seed: random seed for python/numpy    #
#   6) ...: keyword arguments to pass to Python  #
#           API in keyword-value pairs           #
#                                                #
# Output:                                        #
#   1) negative log-likelihood value             #
##################################################
neg_log_lik <- function(pars, parnames, ev_fun, lik_fun, dat, type, ffreq) {
  # compute log-likelihood
  ll <- sum(lik_fun(pars, parnames, ev_fun, dat, type, ffreq))
  -1*ll # return negative of log-likelihood
}

##################################################
# log_post: compute log-posterior value          #
#           for plausible values, return -Inf    #
#           for values ruled out by priors       #
#                                                #
# Inputs:                                        #
#   1) pars: vector of parameter values          #
#   2) parnames: vector of parameter names       #
#   3) priors: list with parnames as names       #
#              and each entry is a list with     #
#              values FUN (density function)     #
#              and param1 and param2 values      #
#   4) model: model function object              #
#   5) dat: vector of observation counts         #
#   6) args: named list of required positional   #
#            arguments, which must include:      #
#       a) n_col: number of columns in grid      #
#       b) n_row: number of rows in grid         #
#       c) mem_length: length of agent memory    #
#       d) seed: random seed for python/numpy    #
#   7) ...: keyword arguments to pass to Python  #
#           API in keyword-value pairs           #
#                                                #
# Output:                                        #
#   1) log-posterior value                       #
##################################################
log_post <- function(pars, parnames, priors, ev_fun, lik_fun, dat, type, ffreq) {
  lp <- log_prior(pars, parnames, priors) # compute log-prior value
  # if the log-prior is -Inf, don't bother computing log-likelihood and return -Inf
  if (!is.finite(lp)) return(-Inf)
  # compute log-likelihood
  ll <- sum(lik_fun(pars, parnames, ev_fun, dat, type, ffreq))
  ll + lp # return log-posterior value
}
