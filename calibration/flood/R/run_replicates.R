library(reticulate)
library(pryr)
library(reshape2)
library(ggplot2)
library(extRemes)
library(scales)

source('R/utils.R')

n_runs <- 1000

river_params <- list('loc' = 865, 'shape' = 0.02, 'scale' = 11)

data_path <- 'data'
dat <- readRDS(file.path(data_path, 'flood_data.rds'))

elev <- dat$elev

river_hist <- dat$hist

models <- c('simple', 'complex')

parnames <- vector('list', length(models))
parnames[['simple']] <- c('int', 'nghd_coef', 'fill_prob')
parnames[['complex']] <- c(parnames[['simple']], 'self_coef')

cols <- 10
rows <- 10
yrs <- 50

start_year <- 2018-yrs+1
end_year <- 2018
mem_length <- 10

out <- vector('list', length(models))
names(out) <- models

for (mod in models) {
  int_true <- -6  # base survival probability
  scoef_true <- 20   # probability of scouting
  ncoef_true <- 4
  fill_true <- 0.01   # probability of surviving a scouting trip
  
  model_params <- list('int'=int_true,
                     'self_coef'=scoef_true,
                     'nghd_coef'=ncoef_true,
                     'fill_prob'=fill_true)

  # add examples and scripts directories to search path
  py_run_string("import sys ")
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

  # source model file to obtain realization
  source_python(file.path(dirname(dirname(getwd())),
                'examples','flood', paste0('flood_', mod, '.py')))
                
  
  model <- partial(BasicFloodModel,
          grid_sz=tuple(as.integer(cols), as.integer(rows)),
          mem_length=as.integer(mem_length),
          river_hist=river_hist,
          elev=elev,
          query=dict('states'=get_states),
          start_year=as.integer(start_year),
          end_year=as.integer(end_year),
          move_params=model_params
       )
       
  out[[mod]] <- vector('list', n_runs)
  
  for (i in 1:n_runs) {
    # evaluate model
    model_out <- run_model(model)(seed=as.integer(i))
    out[[mod]][[i]] <- aperm(model_out$states, c(2, 3, 1))
  }
}

# save replications
saveRDS(out, 'model_runs.rds')

##################### plot spaghetti plot ###################################
# melt output for spaghetti plot
out_vac <- vector('list', length(models))
names(out_vac) <- models
for (mod in models) {
  out_vac[[mod]] <- lapply(out[[mod]][1:100], function(l) {
      data.frame(vac=((cols*rows) - apply(l, 3, sum)),
                 yr=1:dim(l)[3])
  })
}

out_melt <- melt(out_vac, id.vars='yr', measure.vars='vac')
colnames(out_melt)[match('L1', colnames(out_melt))] <- 'model'
colnames(out_melt)[match('L2', colnames(out_melt))] <- 'run'

p <- ggplot(out_melt[out_melt$model == 'complex',]) +
     geom_line(aes(x=yr, y=value, group=run), alpha=0.6, color='grey') +
     geom_line(data=out_melt[out_melt$model == 'complex' & out_melt$run == '1',], aes(x=yr, y=value), alpha=0.6, color='black') +
     scale_x_continuous('Simulation Year', breaks=pretty_breaks()) +
     scale_y_continuous('Number of Vacancies', breaks=pretty_breaks()) +
     theme_bw(base_size=12) +
     theme(panel.grid.major=element_blank(),
           panel.grid.minor=element_blank())
           
pdf('figures/complex_spaghetti.pdf', height=2.5, width=2.5)
p
dev.off()

png('figures/complex_spaghetti.png', height=2.5, width=2.5, units='in', res=300)
p
dev.off()

################ plot model comparison plot #################################
# melt output
out_fin_freq <- lapply(out, function(l) Reduce('+', lapply(l, function(ll) ll[,,yrs]))/length(l))

out_melt <- melt(lapply(out_fin_freq, function(l) data.frame(
              freq=as.numeric(l),
              rl=as.numeric(1/(1-pevd(elev, type='GEV', loc=river_params$loc, shape=river_params$shape, scale=river_params$scale))))),
            id.vars='rl', measure.vars='freq')
colnames(out_melt)[match('L1', colnames(out_melt))] <- 'model'

out_melt$model <- factor(out_melt$model, c('simple', 'complex'))

base_breaks <- function(n = 10){
    function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
}

p <- ggplot(out_melt) +
     geom_line(aes(x=rl, y=value, color=model, linetype=model)) +
     scale_x_log10('Return Period (yr)', breaks=base_breaks(n=5), labels=prettyNum) +
     scale_y_continuous('Occupancy Probability', breaks=pretty_breaks()) +
     theme_bw(base_size=12) +
     theme(panel.grid.major=element_blank(),
           panel.grid.minor=element_blank()) +
     scale_color_brewer('Model', palette='Set2', labels=c('No Interactions', 'Spatial Interactions')) +
     scale_linetype('Model', labels=c('No Interactions', 'Spatial Interactions'))
     
png('figures/model_comparison.png', height=2.5, width=5, units='in', res=300)
p
dev.off()

pdf('figures/model_comparison.pdf', height=2.5, width=5)
p
dev.off()
