# This script plots hindcasts for the flood vacancy agent-based model

library(reticulate) # interface with Python
library(ncdf4)  # interface with netcdf files
library(reshape2)
library(pryr)
library(foreach)
library(doParallel)

source('utils.R')

set.seed(1234)

# set number of prior and posterior samples
nsamp <- 50000

args <- commandArgs(trailingOnly=TRUE)
lik_type <- args[1]

data_path <- '../data'
out_path <- file.path('../output', paste0('model-id-', lik_type))

# if PBS array used, get job number and specify parameters
# otherwise get arguments from command call
aid <- Sys.getenv('PBS_ARRAYID')
if (aid != '') {
  id <- as.numeric(aid)
  n_yrs <- c(10, 25, 50)
  n_rows <- c(5, 10)
  n_cols <- c(5, 10)
  priors <- c('informed', 'weakly')
  models <- c('simple', 'complex')
  
  pcases <- expand.grid(yr=n_yrs, cols=n_cols, rows=n_rows, prior=priors, model=models)
  yrs <- pcases[id, 'yr']
  rows <- pcases[id, 'rows']
  cols <- pcases[id, 'cols']
  prior <- pcases[id, 'prior']
  model <- pcases[id, 'model']
} else {
  yrs <- as.numeric(args[2])
  rows <- as.numeric(args[3])
  cols <- as.numeric(args[4])
  prior <- args[5]
  model <- args[6]
}

start_year <- 2018-yrs+1
end_year <- 2018
mem_length <- 10

# load positive control data
pcout_file <- file.path(data_path, 'model-id',
                          paste('pcout-', yrs, '-',
                                  rows, '-', cols, '.nc',
                                  sep = ''
                               )
                        )
f <- nc_open(pcout_file)
elev <- ncvar_get(f, 'elev')
river_hist <- ncvar_get(f, 'water_height')
nc_close(f)

# loop over cases and estimate prior and posterior predictive distributions
# register cluster
print('Starting cluster...')
# start cluster
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

hind_path <- file.path('../hind', paste0('model-id-', lik_type))
if (!dir.exists(hind_path)) {
  dir.create(hind_path)
}

dists <- c('prior', 'posterior')

hind_out <- foreach(i=1:length(dists),
              .packages=c('reticulate', 'pryr'),
              .export=c('elev', 'river_hist', 'run_model',
                'nsamp', 'dists', 'yrs', 'rows', 'cols', 'prior', 'model',
                'start_year', 'end_year', 'mem_length', 'out_path')
            ) %dopar% {
  # source Python model script
  # add examples and scripts directories to search path
  py_run_string("import sys")
  py_run_string(paste0("sys.path.append('",
                  file.path(dirname(dirname(getwd())), 'examples', 'flood'),
                  "')"
                )
  )
  py_run_string(paste0("sys.path.append('",
                  file.path(dirname(dirname(getwd())), 'scripts'),
                  "')"
               )
  )
                  
  # source model file
  source_python(file.path(dirname(dirname(getwd())),
                'examples','flood', paste0('flood_', model, '.py')))
  
  
  # define partial model run function
  mod <- partial(BasicFloodModel,
              grid_sz=tuple(as.integer(cols), as.integer(rows)),
              mem_length=as.integer(mem_length),
              seed=as.integer(1),
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
    mcmc_file <- file.path(out_path, model, prior,
    	                      paste('MCMC-', yrs, '-',
           			            rows, '-', cols, '.rds',
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
  paste('hind-', yrs, '-', rows, '-', cols, '-', model, '-', prior, '.rds', sep=''))

saveRDS(hind_out, hind_file)

print('Done!')