library(reshape2)
library(plyr)
library(ggplot2)

yrs <- c(10, 25, 50)
cols <- c(5, 10)
rows <- c(5, 10)

cases <- expand.grid(years=yrs, cols=cols, rows=rows)
# remove extra case
cases <- cases[!(cases$cols == 5 & cases$rows == 10),]

lt <- c('ind', 'agg')

bf <- vector('list', 2)
names(bf) <- lt
prior <- 'informed'

compute_bf <- function(yrs, rows, cols, lik_type, prior) {
  ml_path <- file.path('../ml', paste0('model-id-', lik_type)) # set path for ml output
  # load the ml output file and compute the bayes factor
  ml_file <- file.path(ml_path, paste('ml-', yrs, '-', rows, '-', cols, '-', prior, '.rds', sep=''))
  if (!file.exists(ml_file)) {return(NA)}
  ml <- readRDS(ml_file)
  # compute and return bayes factor
  round(ml$complex$ml - ml$simple$ml, 2)

}

for (lik_type in lt) {
  bf[[lik_type]] <- cases
  bf[[lik_type]]$bf <- mapply(compute_bf, cases$years, cases$rows, cases$cols, lik_type=lik_type, prior=prior)
}

bf <- lapply(bf, function(x) cbind(x, agents=x$rows*x$cols))
bf.melt <- melt(bf, id=c('years', 'rows', 'cols', 'agents'))
bf.melt$agents <- as.factor(bf.melt$agents)
#bf.melt$years <- as.factor(bf.melt$years)
colnames(bf.melt)[length(colnames(bf.melt))] <- 'lt'

bf.levels <- data.frame(val=c(0, 1, 3, 5), lab=c('none', 'positive', 'strong', 'very strong'))

lt_labels <- c('ind'='Individual Parcel', 'agg'='Aggregated')

# plot
p <- ggplot() +
  geom_point(data=bf.melt[bf.melt$agents != '25',],
    aes(x=years, y=value, color=agents, shape=agents), size=2.5) +
  geom_hline(data=bf.levels, aes(yintercept=val), linetype='dashed') +
  geom_text(data=bf.levels, aes(y=(val+0.3)*1.25, x=40, label=lab), size=3) +
  scale_y_log10('log-Bayes Factor', limits=c(0.1, 70),
    minor_breaks=c(0.5, 5, 50)) +
  scale_x_continuous('Years') +
  scale_color_brewer('Parcels', palette='Set1') +
  scale_shape('Parcels') +
  theme_bw(base_size=12) +
  theme(legend.position='bottom',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  facet_wrap(vars(lt), nrow=2, labeller=labeller(lt=lt_labels), strip='left')

pdf(paste0('figures/bf.pdf'), height=5, width=3)
p
dev.off()

png(paste0('figures/bf.png'), height=5, width=3, units='in', res=300)
p
dev.off()

#png(file.path('figures', 'poster', 'bf.png'), height=6.5, width=13, units='in', res=300)
#p
#dev.off()
th=13, units='in', res=300)
#p
#dev.off()
