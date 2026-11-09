# One final ES figure per estimable trait, matching the pressure-panel reference.
# No custom or anonymous functions. Run from the repository root.
library(ggplot2)
stopifnot(file.exists('reference/ES_verified.csv'))
es <- read.csv('reference/ES_verified.csv')
es <- es[es$Level=='Library73'&es$Weighting=='CF3_matched',]
meta <- read.csv('data/sample_env.csv',encoding='UTF-8')
dictionary <- read.csv('data/trait_dictionary.csv',encoding='UTF-8')
ix <- match(es$Unit_ID,meta$sample_id)
stopifnot(nrow(es)==73*42,!anyNA(ix))
es$Regime <- factor(meta$regime[ix],levels=c('Intermittent','Perennial'))
es$Pressure <- factor(meta$pressure[ix],levels=c('Low','High'),labels=c('Low pressure','High pressure'))
out <- 'ES_pressure_panel_figures'
dir.create(out,showWarnings=FALSE)
index <- list()
figure <- 0L
for(tr in dictionary$Column) {
 z <- es[es$Trait==tr&is.finite(es$ES),]
 if(nrow(z)==0)next
 stopifnot(nrow(z)==73)
 figure <- figure+1L
 # Both pressure panels share the same symmetric y-scale. Most traits use -3..3,
 # with major ticks at -2, 0 and 2. Extend only when needed to retain every point.
 largest <- max(abs(z$ES))
 limit <- if(largest<=3)3 else ceiling(largest/2)*2
 breaks <- seq(-2*floor(limit/2),2*floor(limit/2),by=2)
 title <- gsub('_',' ',tr)
 p <- ggplot(z,aes(Regime,ES))+
  geom_hline(yintercept=0,colour='black',linewidth=.65)+
  geom_boxplot(width=.38,outlier.shape=NA,fill='white',colour='grey20',linewidth=.65,alpha=.6)+
  geom_point(position=position_jitter(width=.08,height=0,seed=20261109),
             colour='grey25',size=1.8,alpha=.6)+
  facet_wrap(~Pressure,nrow=1,drop=FALSE)+
  scale_y_continuous(limits=c(-limit,limit),breaks=breaks,expand=expansion(mult=.025))+
  scale_x_discrete(drop=FALSE)+
  labs(title=title,x='Hydrological regime',y='Standardized effect size (ES)')+
  theme_grey(base_size=14)+theme(plot.title=element_text(hjust=.5,size=18),
   strip.text=element_text(size=13),panel.grid.minor.x=element_blank(),
   axis.text.x=element_text(size=12),plot.margin=margin(10,12,10,12))
 stem <- sprintf('ES_%02d_%s_PRESSURE_PANELS',figure,tr)
 ggsave(file.path(out,paste0(stem,'.png')),p,width=8.5,height=6.5,dpi=400,bg='white')
 ggsave(file.path(out,paste0(stem,'.tiff')),p,width=8.5,height=6.5,dpi=600,compression='lzw',bg='white')
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=8.5,height=6.5,bg='white')
 index[[figure]] <- data.frame(Figure=figure,Trait=tr,Title=title,Filename=stem,
  Profiles=nrow(z),Lower_limit=-limit,Upper_limit=limit,Observed_min=min(z$ES),Observed_max=max(z$ES))
 cat(figure,tr,'range',-limit,'to',limit,'\n');flush.console()
}
manifest <- do.call(rbind,index)
stopifnot(nrow(manifest)==40,all(manifest$Observed_min>=manifest$Lower_limit),
 all(manifest$Observed_max<=manifest$Upper_limit))
write.csv(manifest,file.path(out,'figure_index.csv'),row.names=FALSE)
writeLines(c('Stellate_colony: ES not estimable because null variance is zero.',
 'Other_colony: ES not estimable because null variance is zero.'),file.path(out,'undefined_categories.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
