library(tmvtnorm)
library(batchmeans)

# starting value is not quite the reciprocal importance sampling estimate
# from Gelfand and Dey (1994)
# due to infinite values obtained when exponentiating, we start with the reciprocal mean of
recip.imp.samp <- function(log.p,log.imp) {
  log.ratio <- log.imp - log.p
  -mean(log.ratio)
}

# function to update the bridge sampling estimate at each iteration
# norm.const is the log normalizing constant estimate from the previous iteration
# post, imp are lists with log.p and log.imp passing the associated log-likelihoods
bridge.samp.iter <- function(log.norm.const,
                        post,
                        imp) {

  # normalize posterior likelihoods based on previous normalizing constant estimate
  # some samples (mostly importance samples) might have infinite posterior log-likelihoods
  post.log.p.norm <- post$log.p[is.finite(post$log.p)] - log.norm.const
  imp.log.p.norm <- imp$log.p[is.finite(imp$log.p)] - log.norm.const

  post.log.imp <- post$log.imp[is.finite(post$log.p)]
  imp.log.imp <- imp$log.imp[is.finite(imp$log.p)]

  # get number of samples
  post.num <- length(post.log.p.norm)
  imp.num <- length(imp.log.p.norm)

  # compute updated estimate numerator and denominator
  imp.mean <- mean(exp(imp.log.p.norm)/(imp.num*exp(imp.log.imp)+post.num*exp(imp.log.p.norm)))
  post.mean <- mean(exp(post.log.imp)/(imp.num*exp(post.log.imp)+post.num*exp(post.log.p.norm)))

  # return updated estimate
  log.norm.const + log(imp.mean) - log(post.mean)
}

bridge.samp.rel.err <- function(log.norm.const,
                                post,
                                imp) {

  # normalize posterior likelihoods based on previous normalizing constant estimate
  # some samples (mostly importance samples) might have infinite posterior log-likelihoods
  post.log.p.norm <- post$log.p[is.finite(post$log.p)] - log.norm.const
  imp.log.p.norm <- imp$log.p[is.finite(imp$log.p)] - log.norm.const

  post.log.imp <- post$log.imp[is.finite(post$log.p)]
  imp.log.imp <- imp$log.imp[is.finite(imp$log.p)]

  # get number of samples
  post.num <- length(post.log.p.norm)
  imp.num <- length(imp.log.p.norm)
  imp.factor <- imp.num/(imp.num+post.num)
  post.factor <- post.num/(imp.num+post.num)

  # compute bridging function process estimates for both sequences
  post.bridge.int <- exp(post.log.imp)/(imp.factor*exp(post.log.imp)+post.factor*(exp(post.log.p.norm)))
  imp.bridge.int <- exp(imp.log.p.norm)/(imp.factor*exp(imp.log.imp)+post.factor*(exp(imp.log.p.norm)))

  # return squared relative error estimate
  (var(imp.bridge.int)/mean(imp.bridge.int)^2)/imp.num + (var(post.bridge.int)/mean(post.bridge.int)^2)/post.num
}

# compute bridge sampling estimate of the marginal likelihood relative to a normal importance density
bridge_sample <- function(mcmc_out, tol, npost, nimp, lik_fun, ...) {
  
  # draw posterior samples and get posterior log-likelihood from MCMC output
  post <- list()
  samp_idx <- sample(1:nrow(mcmc_out$samples), npost, replace=T)
  post$samples <- mcmc_out$samples[samp_idx,]
  post$log.p <- mcmc_out$log.p[samp_idx]
  
  # compute importance density mean and covariance matrix using MCMC output
  imp_mean <- bmmat(mcmc_out$samples)[,1]
  imp_cov <- cov(mcmc_out$samples)
  # draw importance samples
  imp <- list()
  lower=c(rep(-Inf, ncol(mcmc_out$samples)-1), 0)
  upper=c(rep(Inf, ncol(mcmc_out$samples)-1), 1)
  imp$samples <- rtmvnorm(n=nimp, lower=lower, upper=upper, mean=imp_mean, sigma=imp_cov)
  
  # compute importance log-likelihood of importance and posterior samples
  imp$log.imp <- dtmvnorm(imp$samples, lower=lower, upper=upper, mean=imp_mean, sigma=imp_cov, log=T)
  post$log.imp <- dtmvnorm(post$samples, lower=lower, upper=upper, mean=imp_mean, sigma=imp_cov, log=T)
  
  # compute posterior log-likelihood of importance samples
  imp$log.p <- apply(apply(imp$samples, 1, lik_fun, ...), 2, sum)
  
  # initialize storage for estimates
  ml <- mat.or.vec(nr=1000,nc=1)
  
  # initialize iterator using reciprocal importance sampling
  ml[1] <- recip.imp.samp(post$log.p, post$log.imp)
  ml[2] <- bridge.samp.iter(ml[1], post[c('log.p','log.imp')],
    imp[c('log.p','log.imp')])
    # iterate until within tolerance.
  t <- 2
  while (abs(ml[t] - ml[t-1]) >= tol) {
    ml[t+1] <- bridge.samp.iter(ml[t], post[c('log.p', 'log.imp')],
      imp[c('log.p', 'log.imp')])
    t <- t+1
  }
  
  ml <- ml[1:t]
  
  # compute the relative standard error of the bridge sampling estimator
  # we can treat the posterior samples as iid due to re-sampling from the posterior,
  # so we use the error formula from Fruhwirth-Schnatter (2004) with the spectral density
  # at frequency 0 set equal to 1.

  re.sq <- bridge.samp.rel.err(ml[length(ml)], post[c('log.p','log.imp')], imp[c('log.p','log.imp')])

  # return ml and re
  list(ml=ml[length(ml)], se=re.sq)
}