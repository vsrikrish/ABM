# This script plots hindcasts for the flood vacancy agent-based model
library(parallel)
library(foreach)
library(doParallel)


source('R/utils.R')
source('R/likelihood.R')

set.seed(1234)

# set number of prior and posterior samples
nsamp <- 5000

# if PBS array used, get job number and specify parameters
# otherwise get arguments from command call
aid <- Sys.getenv('PBS_ARRAYID')
if (aid != '') {
  id <- as.numeric(aid)
  n_yrs <- c(10, 25, 50)
  n_rows <- c(5, 10)
  n_cols <- c(5, 10)
  priors <- c('informed', 'weakly')
  lik_types <- c('agg', 'ind')
  
  pcases <- expand.grid(yr=n_yrs, cols=n_cols, rows=n_rows, prior=priors, lt=lik_types)
  yrs <- pcases[id, 'yr']
  rows <- pcases[id, 'rows']
  cols <- pcases[id, 'cols']
  prior <- pcases[id, 'prior']
  lik_type <- pcases[id, 'lt']
} else {
  args <- commandArgs(trailingOnly=TRUE)
  yrs <- as.numeric(args[1])
  rows <- as.numeric(args[2])
  cols <- as.numeric(args[3])
  prior <- args[4]
}


out_path <- file.path(getwd(), 'output', paste0('model-id-', lik_type))

models <- c('simple', 'complex')

# set up parallel cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

if (lik_type == 'ind') {
  n <- yrs*cols*rows
} else {
  n <- yrs
}

lik_fun <- match.fun(paste0('log_lik_', lik_type))

ll <- foreach(i=1:length(models),
              .packages=c('plyr'),
              .export=c('models', 'ev_mod', 'inv_logit', 'lik_fun', 'yrs', 'rows', 'cols', 'n')
      ) %dopar% {
      
  # load MCMC output

  mcmc_file <- file.path(out_path, models[i], prior,
    	                      paste('MCMC-', yrs, '-',
           			            rows, '-', cols, '.rds',
	                          sep = ''
                          )
                  )
  mcmc <- readRDS(mcmc_file)
  
  # sample from posterior distribution

  idx <- sample(1:nrow(mcmc$mcmc_out$samples), nsamp, replace=TRUE)
  
  # compute log-likelihoods
  vapply(idx, function(j) as.numeric(lik_fun(mcmc$mcmc_out$samples[j,],
            parnames=mcmc$parnames,
            ev_fun=ev_mod,
            dat=mcmc$obs,
            type=models[i],
            ffreq=mcmc$ffreq)),
          numeric(n)
  )
}
names(ll) <- models

# compute expected log-predictive density for each point
lppd <- lapply(ll, function(v) log(rowMeans(exp(v))))
p <- lapply(ll, function(v) apply(v, 1, function(d) var(d)))
elpd <- lapply(models, function(mod) {-2*(lppd[[mod]] - p[[mod]])})
names(elpd) <- models

# estimate waic, the difference, and the standard error of the difference
waic <- lapply(elpd, sum)
waic$delta <- waic[['simple']] - waic[['complex']]
waic$delta_se <- sqrt(n*var(elpd[['simple']]-elpd[['complex']]))

stopCluster(cl)


print('Saving...')
# set location and filename for bridge sampling output
waic_path <- file.path('waic', paste0('model-id-', lik_type))
if (!dir.exists(waic_path)) {
   dir.create(waic_path)
}
waic_file <- file.path(waic_path,
  paste('waic-', yrs, '-', rows, '-', cols, '-', prior, '.rds', sep=''))

saveRDS(waic, waic_file)

print('Done!')
