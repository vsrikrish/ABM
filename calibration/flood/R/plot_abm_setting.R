library(reshape2)
library(ggplot2)
library(extRemes)
library(RColorBrewer)
library(ggpubr)

river_params <- list('loc' = 865, 'shape' = 0.02, 'scale' = 11)

data_path <- 'data'
dat <- readRDS(file.path(data_path, 'flood_data.rds'))

d_rp <- melt(dat$elev)
colnames(d_rp) <- c('col', 'row', 'elevation')
d_rp$return <- 1/(1-pevd(d_rp$elevation, type='GEV', loc=river_params$loc, shape=river_params$shape, scale=river_params$scale))

cols <- brewer.pal(11, 'RdBu')
cols.rect <- brewer.pal(3, 'Accent')

p1 <- ggplot(d_rp, aes(xmin=(col-1), xmax=col, ymin=(row-1), ymax=row, fill=return)) +
geom_rect(col='black') +
scale_fill_gradientn('Return Period \n (yrs)', colors=cols, trans='log', breaks=c(20, 50, 100, 1000)) + theme_bw(base_size=12) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position='bottom', legend.text=element_text(angle=45), legend.key.width=unit(0.25, 'in')) + geom_rect(xmin=0, xmax=10, ymin=0, ymax=10, col=cols.rect[1], fill=alpha('grey', 0), linetype='dashed') + geom_rect(xmin=0, xmax=10, ymin=0, ymax=5, col=cols.rect[2], fill=alpha('grey', 0), linetype='dashed') + geom_rect(xmin=0, xmax=5, ymin=0, ymax=5, col=cols.rect[3], fill=alpha('grey', 0), linetype='dashed')

d_h <- melt(dat$hist)

rl <- list('100'=qevd(0.99, type='GEV', loc=river_params$loc, shape=river_params$shape, scale=river_params$scale), '500'=qevd(0.998, type='GEV', loc=river_params$loc, shape=river_params$shape, scale=river_params$scale), '50'=qevd(0.98, type='GEV', loc=river_params$loc, shape=river_params$shape, scale=river_params$scale))

col.rl <- brewer.pal(4, 'RdBu')

p2 <- ggplot() + geom_line(data=data.frame(value=d_h[1:51,], year=1:51), aes(y=value, x=year), alpha=0.5) + geom_line(data=data.frame(value=d_h[51:100,], year=51:100), aes(y=value, x=year)) + scale_y_continuous('River Maximum Height \n (ft above sea level)') + scale_x_continuous('Simulation Year', breaks=seq(50, 100, by=10), labels=seq(0, 50, by=10)) + theme_bw(base_size=12) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line()) + geom_hline(yintercept=rl$'100', col=col.rl[2], linetype='dashed') + geom_hline(yintercept=rl$'50', col=col.rl[1], linetype='dashed') + geom_hline(yintercept=rl$'500', col=col.rl[4], linetype='dashed') +
geom_text(aes(x=30, y=(rl$'100'+3)), label='100-year return level', size=3) +
geom_text(aes(x=30, y=(rl$'50'+3)), label='50-year return level', size=3) +
geom_text(aes(x=30, y=(rl$'500'+3)), label='500-year return level', size=3)

#p1 <- arrangeGrob(p1, top=textGrob('a', x=unit(0.01, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=24, fontfamily='sans', fontface='bold')))
#p2 <- arrangeGrob(p2, top=textGrob('b', x=unit(0.01, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=24, fontfamily='sans', fontface='bold')))

q <- ggarrange(plotlist=list(p1, p2), ncol=2, nrow=1, widths=c(1, 1), labels='auto', font.label=list(size=12, color='black', face='bold'), legend='bottom', common.legend=FALSE)

#pdf('setting.pdf', height=3.5, width=6)
#grid.arrange(p1, p2, widths=c(1, 1), nrow = 1)
#dev.off()
png(file.path('figures', 'setting.png'), height=3, width=6, res=300, units='in')
print(q)
dev.off()

pdf(file.path('figures', 'setting.pdf'), height=3, width=6)
print(q)
dev.off()
ff()
ff()
