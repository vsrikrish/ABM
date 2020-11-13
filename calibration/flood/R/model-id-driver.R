## This model calibrates the bird ecology model from Thiele et al (2014) using MCMC on a model for the population growth rate

# load libraries
library(reticulate) # interface with Python
library(adaptMCMC)
library(extRemes)
library(doParallel)
library(DEoptim)

use_condaenv("mcmc", required=TRUE)

# load utility functions
source('R/utils.R')
# load density functions
source('R/likelihood.R')

args <- commandArgs(trailingOnly=TRUE)
# if PBS array used, get job number and specify parameters
# otherwise get arguments from command call
aid <- Sys.getenv('PBS_ARRAYID')
if (aid != '') {
  id <- as.numeric(aid)
  n_yrs <- c(10, 25, 50)
  n_rows <- c(5, 10)
  n_cols <- c(5, 10)
  lik <- c('ind', 'agg')
  
  pcases <- expand.grid(yr=n_yrs, cols=n_cols, rows=n_rows, lik=lik, seed=1:5)
  yrs <- pcases[id, 'yr']
  rows <- pcases[id, 'rows']
  cols <- pcases[id, 'cols']
  lik_type <- pcases[id, 'lik']
  seed <- pcases[id, 'seed']
} else {
  yrs <- as.numeric(args[2])
  rows <- as.numeric(args[3])
  cols <- as.numeric(args[4])
  seed <- as.numeric(args[5])
}
# set seed for same background information generation
set.seed(1234)

# set parameters for both data generation and MCMC
start_year <- 2018-yrs+1
end_year <- 2018
mem_length <- 10

# set model names
model_names <- c('simple', 'complex')

data_path <- file.path(getwd(), 'data')
# run positive control test if necessary
if (!dir.exists(data_path)) {
  dir.create(data_path)
}

print('Generating positive control data...')

# run positive control test if it hasn't already been run
pcout_file <- file.path(data_path, 'model-id',
                          paste('pcout-', yrs, '-',
                                  rows, '-', cols, '-', seed, '.rds',
                                  sep = ''
                               )
                        )
# generate water height history and housing elevations
slope <- list('offset' = 890, 'mean' = 10)
elev <- t(slope$offset +
            slope$mean * matrix(seq(0, rows - 1),
                                ncol = cols,
                                nrow = rows
            ) +
            rnorm(cols * rows, 0, 5)
)

river_params <- list('loc' = 865, 'shape' = 0.02, 'scale' = 11)
river_hist <- revd(end_year-start_year+51, loc=river_params$loc, shape=river_params$shape, scale=river_params$scale, type='GEV')
# set "true" parameter values
int_true <- -6  # base survival probability
scoef_true <- 20   # probability of scouting
ncoef_true <- 4
fill_true <- 0.01   # probability of surviving a scouting trip

if (!file.exists(pcout_file)) {

  # add examples and scripts directories to search path
  pysys <- import('sys')
  functools <- import('functools')
  
  pysys$path <- c(pysys$path, file.path(dirname(dirname(getwd())), 'examples', 'flood'), 
                  file.path(dirname(dirname(getwd())), 'scripts'))

  # source complex model file to generate positive control data
  source_python(file.path(dirname(dirname(getwd())),
                'examples','flood', 'flood_complex.py'))

  # add partial inputs to model function
  model <- functools$partial(BasicFloodModel,
            grid_sz=tuple(as.integer(cols), as.integer(rows)),
            mem_length=as.integer(mem_length),
            seed=as.integer(seed),
            river_hist=river_hist,
            elev=elev,
            query=dict('vac'=get_vac,
                       'states'=get_states,
                       'flood_freq'=get_floods),
            start_year=as.integer(start_year),
            end_year=as.integer(end_year)
         )


  # evaluate model
  model_params <- list('int'=int_true,
                       'self_coef'=scoef_true,
                       'nghd_coef'=ncoef_true,
                       'fill_prob'=fill_true)

  model_out <- run_model(model)(move_params=model_params)
  model_out$flood_freq <- aperm(model_out$flood_freq, c(2, 3, 1))
  model_out$states <- aperm(model_out$states, c(2, 3, 1))

  # save to file
  saveRDS(model_out, pcout_file)

    # save vacancy observations
  if (lik_type == 'ind') {
    obs <- model_out$states
  } else if (lik_type == 'agg') {
    obs <- model_out$vac
  }
} else {
  # open pseudoobservation file and pull observation variables
  f <- readRDS(pcout_file)
  if (lik_type == 'ind') {
    obs <- f$states
  } else if (lik_type == 'agg') {
    obs <- f$vac
  }
}

# determine rolling number of floods
floods <- aperm(apply(elev, c(1, 2), function(h) h < river_hist), c(2, 3, 1))
ffreq <- aperm(
    apply(floods, c(1, 2), rollapply, n=mem_length, f=sum),
    c(2,3,1)
  )
ffreq <- ffreq[,,(mem_length / 2 + 1):dim(ffreq)[3]]
ffreq <- ffreq[,,tail(seq_along(1:dim(ffreq)[3]), yrs)]

print('Setting up cluster...')
# set up parallel cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

# set parameter name vector list
parnames <- vector('list', 2)
names(parnames) <- model_names
parnames[['simple']] <- c('int', 'self_coef', 'fill')
parnames[['complex']] <- c('int', 'self_coef', 'nghd_coef', 'fill')

# select likelihood function based on type
lik_fun <- match.fun(paste0('log_lik_', lik_type))
# find MLE estimate to start ------------------------------------------------
print('Finding MLE...')
# parallelize over model cases
init_val <- foreach(i=1:length(model_names),
        .packages=c('DEoptim', 'reticulate', 'plyr'),
        .export=c('obs', 'model_names', 'ffreq', 'ev_mod',
                  'neg_log_lik', 'lik_fun', 'inv_logit',
                  'rollapply')
        ) %dopar% {

  # set up directory for output
  out_path <- file.path(getwd(), 'output', paste0('model-id-', lik_type),
    model_names[i])
  if (!dir.exists(out_path)) {
    dir.create(out_path)
  }

  # set up filename for MLE estimate storage
  init_file <- file.path(out_path,
                  paste('init-', yrs, '-', rows, '-', cols, '-', seed, '.rds',
                  sep='')
               )

  # if the MLE file exists, don't need to rerun
  if (!file.exists(init_file)) {
    # set variable lower and upper bounds
    if (model_names[i] == 'simple') {
      lower_bd <- c(-20, 0, 0)
      upper_bd <- c(20, 40, 1)
    } else {
      lower_bd <- c(-20, 0, 0, 0)
      upper_bd <- c(20, 40, 40, 1)
    }
    init_val <- DEoptim(neg_log_lik,
                  lower = lower_bd,
                  upper = upper_bd,
                  control = DEoptim.control(NP = 200,
                    itermax = 100),
                  parnames = parnames[[model_names[i]]],
                  ev_fun = ev_mod,
                  lik_fun = lik_fun,
                  dat = obs,
                  type= model_names[i],
                  ffreq=ffreq
                )$optim$bestmem
    saveRDS(init_val, init_file)
  } else {
    init_val <- readRDS(init_file)
  }
  init_val
}

# run MCMC -----------------------------------------------------------------
  print('Running preliminary MCMC...')
# set priors
foreach(i=1:length(model_names),
        .packages=c('adaptMCMC', 'reticulate', 'plyr', 'pryr'),
        .export=c('obs', 'model_names', 'ffreq', 'ev_mod',
                  'neg_log_lik', 'lik_fun', 'log_prior', 'log_post',
                  'init_val', 'inv_logit', 'rollapply', 'yrs', 'rows',
                  'cols')
        ) %dopar% {
  priors <- vector('list', length(parnames[[model_names[i]]]))

  names(priors) <- model_names
    for (name in parnames[[model_names[i]]]) {
      if (name == 'int') {
        priors[[name]] <- list(FUN=dnorm,
                               params=list('mean'=-7, 'sd'=1)
                          )
      } else if (name == 'self_coef') {
        priors[[name]] <- list(FUN=dnorm,
                               params=list('mean'=19, 'sd'=2)
                          )
      } else if (name == 'nghd_coef') {
        priors[[name]] <- list(FUN=dnorm,
                               params=list('mean'=5, 'sd'=2)
                          )
      } else if (name == 'fill') {
        priors[[name]] <- list(FUN=dbeta,
                               params=list('shape1'=1, 'shape2'=10)
                          )
      }
    }

  print('Running MCMC...')
  # set MCMC parameters
  niter_prelim <- 30000
  niter_prod <- 150000
  accept_rate_many <- 0.234
  accept_rate_few <- 0.44
  accept_rate <- accept_rate_many + (accept_rate_few-accept_rate_many)/length(parnames)

  # define partial function to simplify MCMC calls
  log_post_mcmc <- partial(log_post, parnames=parnames[[model_names[i]]],
                      priors=priors, ev_fun=ev_mod, lik_fun=lik_fun, dat=obs,
                      type=model_names[i], ffreq=ffreq
                    )

  # set up directory for output
  out_path <- file.path(getwd(), 'output', paste0('model-id-', lik_type),
    model_names[i])
  if (!dir.exists(out_path)) {
    dir.create(out_path)
  }

  # run preliminary MCMC chain to find starting point and jump covariance matrix estimate
  prelim_file <- file.path(out_path,
    paste('prelim_MCMC-', yrs, '-', rows, '-', cols, '-', seed, '.rds', sep='')
  )
  if (!file.exists(prelim_file)) {
    amcmc_prelim <- MCMC(log_post_mcmc, niter_prelim,
                      init_val[[i]],
                      adapt = TRUE, acc.rate = accept_rate, gamma = 0.55,
                      list = TRUE,
                      n.start = max(500, round(0.05*niter_prelim))
                    )
    saveRDS(amcmc_prelim, prelim_file)
  } else {
    amcmc_prelim <- readRDS(prelim_file)
  }

  mcmc_file <- file.path(out_path,
    paste('MCMC-', yrs, '-', rows, '-', cols, '-', seed, '.rds', sep=''))
  cov_jump <- amcmc_prelim$cov.jump
  print('Running MCMC...')
  amcmc_out <- MCMC(log_post_mcmc, niter_prod,
                  amcmc_prelim$samples[nrow(amcmc_prelim$samples), ],
                  scale = cov_jump, list = TRUE, adapt = FALSE
               )
  mcmc_list <- list(mcmc_out=amcmc_out,
                priors=priors,
                obs=obs,
                mle=init_val[[i]],
                parnames=parnames[[model_names[i]]],
                ffreq=ffreq
                )
  saveRDS(mcmc_list, mcmc_file)
}
print('Done!')
