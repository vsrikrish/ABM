library(plyr) # convert lists to data frames
library(reshape2) # melt data frames
library(ggplot2)  # plot functionality
library(ggpubr)
library(scales)

var_lbl <- c(
  int = 'Logit Intercept',
  self_coef = 'Flood Coefficient',
  nghd_coef = 'Vacancy Coefficient',
  fill = 'Vacant Fill Probability'
)

args <- commandArgs(trailingOnly=TRUE)
n_yr <- as.numeric(args[1])
n_row <- as.numeric(args[2])
n_col <- as.numeric(args[3])

lik_type <- c('ind', 'agg')

# set directories
data_path <- file.path(getwd(), 'data')
fig_path <- file.path(getwd(), 'figures')

if (!dir.exists(fig_path)) {
  dir.create(fig_path)
}


parnames <- c('int', 'self_coef', 'nghd_coef', 'fill')
vals <- c(-6, 20, 4, 0.01)

priors <- c('informed', 'weakly')

pars_true <- data.frame(variable=parnames, vals=vals)
pars_true$var_f <- factor(pars_true$variable, levels=parnames)
pars_true$var <- plyr::revalue(pars_true$var_f, var_lbl)
samps <- vector('list', length(parnames))
names(samps) <- parnames

p <- vector('list', length(lik_type))
names(p) <- lik_type

n_samp <- 100000



for (pri in priors) {

  for (lt in lik_type) {
    out_path <- file.path(getwd(), 'output', paste0('model-id-', lt))

    # load MCMC output
    mcmc_file <- file.path(out_path, 'complex', pri,
                            paste('MCMC-', n_yr, '-',
                                    n_row, '-', n_col, '.rds',
                                    sep = ''
                                  )
                          )
                          
    mcmc_list <- readRDS(mcmc_file)

    for (name in parnames) {
      samps[[name]] <- vector('list', 2)
      names(samps[[name]]) <- c('prior', 'posterior')
      prior <- mcmc_list$priors[[name]]
      if (name == 'fill') {
        samps[[name]][['prior']] <- do.call(rbeta,
                                      c(list(n=n_samp), prior$params)
                                    )
      } else {
        samps[[name]][['prior']] <- do.call(rnorm,
                                      c(list(n=n_samp), prior$params)
                                    )
      }
      post_idx <- sample(1:nrow(mcmc_list$mcmc_out$samples), n_samp, replace=T)
      samps[[name]][['posterior']]<- mcmc_list$mcmc_out$samples[post_idx,
                                       match(name, parnames)]
    }

    # melt samples
    samps_melt <- melt(samps)
    colnames(samps_melt) <- c('sample', 'distribution', 'variable')
    samps_melt$var_f <- factor(samps_melt$variable, levels=parnames)
    samps_melt$var <- plyr::revalue(samps_melt$var_f, var_lbl)
    samps_melt$distribution <- factor(samps_melt$distribution, c('prior', 'posterior'))

  # plot
    p[[lt]] <- ggplot(samps_melt) +
          stat_density(aes(sample, fill=distribution, group=distribution), geom='area', position='identity', alpha=0.5) +
          geom_vline(data=pars_true, aes(xintercept=vals), color='black', alpha=0.5, linetype='dashed', size=1) +
          facet_wrap(vars(var), ncol=1, scales='free',
            labeller=labeller(var=label_wrap_gen(15)),
            strip = 'right'
          ) +
          scale_fill_brewer('', palette='Dark2', labels=c('Prior', 'Posterior')) +
          scale_x_continuous('Value', breaks=pretty_breaks()) +
          scale_y_continuous('Density', breaks=pretty_breaks()) +
          theme_bw(base_size=12) +
          theme(
            panel.grid.major=element_blank(),   panel.grid.minor=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.spacing.y=unit(0.15, 'lines')
          )
  }
  
  q <- ggarrange(plotlist=p, ncol=2, nrow=1, widths=c(1, 1), labels='auto', font.label=list(size=12, color='black', face='bold'), legend='right', common.legend=TRUE)

  png(file.path(fig_path, paste(pri, '-dist', '.png', sep='')), width=7, height=4, units='in', res=300)
  print(q)
  dev.off()


 pdf(file.path(fig_path, paste(pri, '-dist', '.pdf', sep='')), width=7, height=4)
  print(q)
  dev.off()

}
}