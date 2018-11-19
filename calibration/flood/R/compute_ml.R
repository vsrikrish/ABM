source('R/utils.R')
source('R/likelihood.R')
source('bridge_sample.R')
library(foreach)
library(doParallel)

set.seed(1234)

# set number of posterior and importance samples for the bridge sampling estimate
nsamp <- 5000
# set tolerance for bridge sampling iterator
TOL <- 1e-3

args <- commandArgs(trailingOnly=TRUE)
lik_type <- args[1]
# if PBS array used, get job number and specify parameters
# otherwise get arguments from command call
aid <- Sys.getenv('PBS_ARRAYID')
if (aid != '') {
  id <- as.numeric(aid)
  n_yrs <- c(10, 25, 50)
  n_rows <- c(5, 10)
  n_cols <- c(5, 10)
  priors <- c('informed', 'weakly')
  
  pcases <- expand.grid(yr=n_yrs, cols=n_cols, rows=n_rows, prior=priors)
  yrs <- pcases[id, 'yr']
  rows <- pcases[id, 'rows']
  cols <- pcases[id, 'cols']
  prior <- pcases[id, 'prior']
} else {
  yrs <- as.numeric(args[2])
  rows <- as.numeric(args[3])
  cols <- as.numeric(args[4])
  prior <- args[5]
}

model_names <- c('simple', 'complex')

print('Starting cluster...')
# start cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

if (!dir.exists('ml')) {
  dir.create('ml')
}
ml_path <- file.path('ml', paste0('model-id-', lik_type))
if (!dir.exists(ml_path)) {
  dir.create(ml_path)
}

lik_fun <- match.fun(paste0('log_lik_', lik_type))

print('Estimating marginal likelihoods...')
bs_out <- foreach(i=1:length(model_names),
            .packages=c('tmvtnorm', 'batchmeans', 'plyr'),
            .export=c('yrs','rows', 'cols', 'lik_type', 'prior', 'model_names',
              'ml_path', 'TOL', 'nsamp', 'lik_fun', 'ev_mod',
              'bridge_sample', 'recip.imp.samp', 'bridge.samp.iter',
              'bridge.samp.rel.err', 'inv_logit')
          ) %dopar% {
  # set location and filename for MCMC output
  out_path <- file.path('../output', paste0('model-id-', lik_type), model_names[i], prior)
  mcmc_file <- file.path(out_path,
    paste('MCMC-', yrs, '-', rows, '-', cols, '.rds', sep=''))
  # read in MCMC output
  mcmc <- readRDS(mcmc_file)
  
  # run the bridge sampling
  bridge_sample(mcmc_out=mcmc$mcmc_out, tol=TOL, npost=nsamp, nimp=nsamp, lik_fun=lik_fun, parnames=mcmc$parnames, ev_fun=ev_mod, dat=mcmc$obs, type=model_names[i], ffreq=mcmc$ffreq)
}

stopCluster(cl)

names(bs_out) <- model_names


print('Saving...')
# set location and filename for bridge sampling output
ml_file <- file.path(ml_path,
  paste('ml-', yrs, '-', rows, '-', cols, '-', prior, '.rds', sep=''))

saveRDS(bs_out, ml_file)

print('Done!')
