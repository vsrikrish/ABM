# This script plots hindcasts for the flood vacancy agent-based model

library(reticulate) # interface with Python
library(extRemes)
library(reshape2)
library(pryr)
library(foreach)
library(doParallel)

source('R/utils.R')

use_condaenv("mcmc", conda='~/work/miniconda3/bin/conda', required = TRUE)

set.seed(1234)

# set number of prior and posterior samples
nsamp <- 50000

args <- commandArgs(trailingOnly=TRUE)

# if PBS array used, get job number and specify parameters
# otherwise get arguments from command call
aid <- Sys.getenv('PBS_ARRAYID')
if (aid != '') {
  id <- as.numeric(aid)
  n_yrs <- c(10, 25, 50)
  n_rows <- c(5, 10)
  n_cols <- c(5, 10)
  models <- c('simple', 'complex')
  liks <- c('ind', 'agg')
  
  pcases <- expand.grid(yr=n_yrs, cols=n_cols, rows=n_rows, model=models, lik=liks, seed=1:5)
  yrs <- pcases[id, 'yr']
  rows <- pcases[id, 'rows']
  cols <- pcases[id, 'cols']
  model <- pcases[id, 'model']
  seed <- pcases[id, 'seed']
  lik_type <- pcases[id, 'lik']
  print(pcases[id,])
} else {
  yrs <- as.numeric(args[1])
  rows <- as.numeric(args[2])
  cols <- as.numeric(args[3])
  model <- args[4]
  seed <- as.numeric(args[5])
  lik_type <- args[6]
}

data_path <- 'data'
out_path <- file.path('output', paste0('model-id-', lik_type))

start_year <- 2018-yrs+1
end_year <- 2018
mem_length <- 10

# load positive control data
pcout_file <- file.path(data_path, 'model-id',
                          paste('pcout-', yrs, '-',
                                  rows, '-', cols, '-', seed, '.rds',
                                  sep = ''
                               )
                        )
f <- readRDS(pcout_file)$flood_freq
# generate water height history and housing elevations
slope <- list('offset' = 890, 'mean' = 10)
elev <- t(slope$offset +
            slope$mean * matrix(seq(0, rows - 1),
                                ncol = cols,
                                nrow = rows
            ) +
            rnorm(cols * rows, 0, 5))
river_params <- list('loc' = 865, 'shape' = 0.02, 'scale' = 11)
river_hist <- revd(end_year-start_year+51, loc=river_params$loc, shape=river_params$shape, scale=river_params$scale, type='GEV')

# loop over cases and estimate prior and posterior predictive distributions
# register cluster
print('Starting cluster...')
# start cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

hind_path <- file.path('hind', paste0('model-id-', lik_type))
if (!dir.exists(hind_path)) {
  dir.create(hind_path)
}

dists <- c('prior', 'posterior')

hind_out <- foreach(i=1:length(dists),
              .packages=c('reticulate', 'pryr'),
              .export=c('elev', 'river_hist', 'run_model',
                'nsamp', 'dists', 'yrs', 'rows', 'cols', 'model', 'seed', 
                'start_year', 'end_year', 'mem_length', 'out_path')
            ) %dopar% {
  # source Python model script
  # add examples and scripts directories to search path
  # add examples and scripts directories to search path
  pysys <- import('sys')
  functools <- import('functools')
  
  pysys$path <- c(pysys$path, file.path(dirname(dirname(getwd())), 'examples', 'flood'), 
                  file.path(dirname(dirname(getwd())), 'scripts'))
              
  
  # source model file
  source_python(file.path(dirname(dirname(getwd())),
                'examples','flood', paste0('flood_', model, '.py')))
  
  
  # define partial model run function
  mod <- partial(BasicFloodModel,
              grid_sz=tuple(as.integer(cols), as.integer(rows)),
              mem_length=as.integer(mem_length),
              seed=as.integer(seed),
              river_hist=river_hist,
              elev=elev,
              query=dict('vac'=get_vac),
              start_year=as.integer(start_year),
              end_year=as.integer(end_year)
  )

  # define wrapper for hindcast model runs
  model_wrap <- function(pars, parnames, model) {
    force(model)
  
    int <- pars[match('int', parnames)]
    scoef <- pars[match('self_coef', parnames)]
    ncoef <- pars[match('nghd_coef', parnames)]
    fill <- pars[match('fill', parnames)]
  
    model_params <- list('int'=int,
                       'self_coef'=scoef,
                       'fill_prob'=fill)
    if (!is.na(ncoef)) { model_params$nghd_coef <- ncoef }
    
    model_out <- run_model(model)(move_params=model_params)
    model_out$vac
  }

  # load MCMC output
    mcmc_file <- file.path(out_path, model,
    	                      paste('MCMC-', yrs, '-',
           			            rows, '-', cols, '-', seed, '.rds',
	                          sep = ''
                          )
                  )
  mcmc <- readRDS(mcmc_file)
                      
  parnames <- mcmc$parnames
  
    # if generating the prior predictive, sample
    if (dists[i] == 'prior') {
      samps <- data.frame(matrix(ncol=length(parnames), nrow=nsamp))
      for (name in parnames) {
        # find parameter value and call prior density function
        if (name == 'fill') {
          samps[, match(name, parnames)] <- do.call(rbeta, c(list(n=nsamp), mcmc$priors[[name]]$params))
        } else {
          samps[, match(name, parnames)] <- do.call(rnorm, c(list(n=nsamp), mcmc$priors[[name]]$params))
        }
      }
    # else sample from the MCMC chain
    } else {
    # generate MCMC ensemble
      idx <- sample(1:nrow(mcmc$mcmc_out$samples), nsamp, replace=T)
      samps <- mcmc$mcmc_out$samples[idx,]
    }

  # generate hindcasts
  apply(samps, 1, model_wrap,
      parnames=mcmc$parnames,
      model=mod
  )

}

stopCluster(cl)

names(hind_out) <- dists

print('Saving...')
# set location and filename for bridge sampling output
hind_file <- file.path(hind_path,
  paste('hind-', yrs, '-', rows, '-', cols, '-', model, '-', seed, '.rds', sep=''))

saveRDS(hind_out, hind_file)

print('Done!')