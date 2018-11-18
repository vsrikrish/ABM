# This script plots hindcasts for the flood vacancy agent-based model

library(ggplot2)  # plot functionality
library(ncdf4)  # interface with netcdf files
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(scales)

# if PBS array used, get job number and specify parameters
# otherwise get arguments from command call
n_yrs <- c(10, 25, 50)
n_rows <- c(5, 10)
n_cols <- c(5, 10)
priors <- c('informed', 'weakly')
models <- c('simple', 'complex')
lik_types <- c('agg', 'ind')

pcases <- expand.grid(yr=n_yrs, cols=n_cols, rows=n_rows, prior=priors, model=models, lt=lik_types)
pcases$agents <- pcases$rows*pcases$cols
agents <- unique(pcases$agents)

# set up storage for hindcasts
hind <- vector('list', length(models))
names(hind) <- models

for (model in models) {
  hind[[model]] <- vector('list', length(n_yrs))
  names(hind[[model]]) <- as.character(n_yrs)
  for (yr in as.character(n_yrs)) {
    hind[[model]][[yr]] <- vector('list', length(agents))
    names(hind[[model]][[yr]]) <- as.character(agents)
    for (agent in as.character(agents)) {
      hind[[model]][[yr]][[agent]] <- vector('list', length(lik_types))
      names(hind[[model]][[yr]][[agent]]) <- lik_types
      for (lt in lik_types) {
        hind[[model]][[yr]][[agent]][[lt]] <- vector('list', length(priors))
        names(hind[[model]][[yr]][[agent]][[lt]]) <- priors
      }
    }
  }
}


obs <- vector('list', length(models))
names(obs) <- models

for (model in models) {
  obs[[model]] <- vector('list', length(n_yrs))
  names(obs[[model]]) <- as.character(n_yrs)
  for (yr in as.character(n_yrs)) {
    obs[[model]][[yr]] <- vector('list', length(agents))
    names(obs[[model]][[yr]]) <- as.character(agents)
    for (agent in as.character(agents)) {
      obs[[model]][[yr]][[agent]] <- vector('list', length(lik_types))
      names(obs[[model]][[yr]][[agent]]) <- lik_types
      for (lt in lik_types) {
        obs[[model]][[yr]][[agent]][[lt]] <- vector('list', length(priors))
        names(obs[[model]][[yr]][[agent]][[lt]]) <- priors
      }
    }
  }
}

data_path <- file.path(dirname(getwd()), 'data')

for (i in 1:nrow(pcases)) {
  yrs <- pcases[i, 'yr']
  rows <- pcases[i, 'rows']
  cols <- pcases[i, 'cols']
  prior <- pcases[i, 'prior']
  mod <- pcases[i, 'model']
  agents <- as.character(pcases[i, 'agents'])
  yr <- as.character(yrs)
  lik_type <- pcases[i, 'lt']

  # set directories
  hind_path <- file.path(dirname(getwd()), 'hind', paste0('model-id-', lik_type))

  if (rows <= cols) {

    # load positive control data
    # load hindcasts
    hind_file <- file.path(hind_path,
      paste('hind-', yrs, '-', rows, '-', cols, '-', mod, '-', prior, '.rds', sep=''))

    hc <- readRDS(hind_file)

    # compute quantiles and medians
    hind[[mod]][[yr]][[agents]][[lik_type]][[prior]] <-
      lapply(hc, function(x) as.data.frame(t(apply(x, 1, quantile, probs=c(0.025, 0.5, 0.975)))))

    hind[[mod]][[yr]][[agents]][[lik_type]][[prior]] <-
      lapply(hind[[mod]][[yr]][[agents]][[lik_type]][[prior]], setNames, c('lower', 'median', 'upper'))

    hind[[mod]][[yr]][[agents]][[lik_type]][[prior]] <-
      lapply(hind[[mod]][[yr]][[agents]][[lik_type]][[prior]], function(x) cbind(x, year=1:yrs))

    pcout_file <- file.path(data_path, 'model-id',
                          paste('pcout-', yrs, '-',
                                  rows, '-', cols, '.nc',
                                  sep = ''
                               )
                        )

    f <- nc_open(pcout_file)
    obs[[mod]][[yr]][[agents]][[lik_type]][[prior]] <- data.frame(obs=ncvar_get(f, 'vac'), year=1:yrs)
    nc_close(f)
  }
}

hind.melt <- melt(hind, id.var=c('year', 'lower', 'median', 'upper'))
obs.melt <- melt(obs, id.var=c('year', 'obs'))

colnames(hind.melt) <- c('year', 'lower', 'median', 'upper', 'dist', 'prior', 'lik_type', 'parcels', 'years', 'model')
colnames(obs.melt) <- c('year', 'obs', 'prior', 'lik_type', 'parcels', 'years', 'model')

hind.melt$ymax <- as.numeric(hind.melt$parcels)
hind.melt$parcels <- factor(hind.melt$parcels, levels(factor(hind.melt$parcels))[c(2,3,1)])
hind.melt$years <- factor(hind.melt$years)
hind.melt$model <- factor(hind.melt$model)
hind.melt$prior <- factor(hind.melt$prior)
hind.melt$dist <- factor(hind.melt$dist, c('prior', 'posterior'))

obs.melt$parcels <- factor(obs.melt$parcels, levels(factor(obs.melt$parcels))[c(2,3,1)])
obs.melt$years <- factor(obs.melt$years)
obs.melt$model <- factor(obs.melt$model)
obs.melt$prior <- factor(obs.melt$prior)

colors <- brewer.pal(3, 'Dark2')
colors[3] <- 'black'
names(colors) <- c('prior', 'posterior', 'obs')

gcases <- expand.grid(models=models, lt=lik_types)
p <- vector('list', nrow(gcases))

for (i in 1:nrow(gcases)) {
  lt <- gcases[i, 'lt']
  mod <- gcases[i, 'models']
  
  hind_sub <- hind.melt[(hind.melt$prior == 'informed') & (hind.melt$model == mod) & (hind.melt$lik_type == lt), ]
  
  p[[i]] <- ggplot() +
      geom_ribbon(data=hind_sub,
        aes(x=year, ymin=lower, ymax=upper, fill=dist), alpha=0.2) +
      geom_line(data=hind_sub[hind_sub$dist == 'posterior',], aes(x=year, y=median, color=dist), size=1) +
      geom_point(data=obs.melt, aes(x=year, y=obs, color='obs'), size=0.75, alpha=0.6) +
      geom_blank(data=hind_sub, aes(x=0, y=ymax)) +
      theme_bw(base_size=12) +
      scale_x_continuous('Simulation Year', breaks=pretty_breaks()) +
      scale_y_continuous('Housing Vacancies', breaks=pretty_breaks()) +
      scale_color_manual('', values=colors, labels=c('Observations', 'Posterior Median')) +
      scale_fill_manual('', values=colors[1:2], labels=c('Prior 95% CI', 'Posterior 95% CI')) +
      facet_grid(parcels ~ years, scales='free', labeller=label_both) +
      theme(legend.position='bottom', legend.box='horizontal',
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.spacing.y=unit(-0.5, 'lines'),
        panel.spacing.x=unit(0.4,'lines'),
        panel.spacing.y=unit(0.25, 'lines')
       ) +
      guides(color=guide_legend(override.aes=list(
        linetype=c(0,1), shape=c(16, NA)), order=1
        ), fill=guide_legend(order=2))


}
  

q <- ggarrange(plotlist=p, ncol=2, nrow=2, heights=c(1, 1), widths=c(1, 1), labels='auto', font.label=list(size=12, color='black', face='bold'), legend='bottom', common.legend=TRUE)

png(paste0('hindcast-informed.png'), width=8, height=8, units='in', res=300)
print(q)
dev.off()

pdf(paste0('hindcast-informed.pdf'),
  width=8, height=8)
print(q)
dev.off()
.off()
