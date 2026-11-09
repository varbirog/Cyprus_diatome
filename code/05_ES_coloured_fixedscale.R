# Final trait-group ES figures from the fixed verified publication null.
# No custom/anonymous functions; run from the repository root.
stopifnot(file.exists('reference/ES_verified.csv'))
library(ggplot2)
d <- read.csv('reference/ES_verified.csv',check.names=FALSE)
d <- d[d$Level=='Library73'&d$Weighting=='CF3_matched',]
meta <- read.csv('data/sample_env.csv',check.names=FALSE,encoding='UTF-8')
stopifnot(nrow(d)==73*42,nrow(meta)==73)
ix <- match(d$Unit_ID,meta$sample_id)
stopifnot(!anyNA(ix))
d$Group <- factor(meta$group[ix],levels=c('Perennial_Low','Perennial_High','Intermittent_Low','Intermittent_High'))
families <- list(
 Size=c('S1','S2','S3','S4','S5','Pioneer'),
 Shape=c('LW1','LW2','LW3','LW4','LW5','LW6'),
 Moisture=c('VSA','AOC','AER','NAE','TER'),
 Attachment=c('Not_attached','Attached','Adnate','Pedunculate','Pad','Stalk'),
 Colony=c('Colonial','Non_colonial','Mucous_tubule_colony','Chain_colony','Zigzag_colony','Rosette_colony','Ribbon_colony','Stellate_colony','Arbuscular_colony','Other_colony'),
 Guild=c('High_profile_guild_v2','Low_profile_guild_v2','Motile_guild_v2','Euplanctonic_guild_v2'),
 Habitat_mobility=c('Benthic','Planctonic','Marine','Freshwater','Mobile'))
titles <- c(Size='Cell size and colonisation',Shape='Cell shape',Moisture='Moisture preference',
 Attachment='Attachment type',Colony='Colony type',Guild='Ecological guilds',Habitat_mobility='Life form habitat and mobility')
stopifnot(setequal(unlist(families),unique(d$Trait)))
colours <- c(Perennial_Low='#0072B2',Perennial_High='#56B4E9',Intermittent_Low='#D55E00',Intermittent_High='#E69F00')
shapes <- c(Perennial_Low=16,Perennial_High=17,Intermittent_Low=15,Intermittent_High=3)
labels <- c('Perennial / Low','Perennial / High','Intermittent / Low','Intermittent / High')
out <- 'ES_coloured_fixedscale_figures'
dir.create(out,showWarnings=FALSE)
captions <- character()
for(i in seq_along(families)) {
 family <- names(families)[i]
 traits <- families[[i]]
 panel <- d[d$Trait %in% traits,]
 panel$Trait_label <- factor(gsub('_',' ',panel$Trait),levels=gsub('_',' ',traits))
 usable <- panel[is.finite(panel$ES),]
 undefined <- panel[!is.finite(panel$ES),]
 undefined <- undefined[!duplicated(undefined$Trait),]
 undefined$Group <- factor(rep('Perennial_High',nrow(undefined)),levels=levels(panel$Group))
 axis_rows <- list()
 for(tr in traits) {
  values <- usable$ES[usable$Trait==tr]
  largest <- if(length(values))max(abs(values))else 0
  limit <- if(largest<=3)3 else ceiling(largest/2)*2
  axis_rows[[length(axis_rows)+1L]] <- data.frame(
   Trait_label=factor(gsub('_',' ',tr),levels=levels(panel$Trait_label)),
   Group=factor('Perennial_Low',levels=levels(panel$Group)),Limit=c(-limit,limit))
 }
 axis_limits <- do.call(rbind,axis_rows)
 p <- ggplot(usable,aes(Group,ES))+
  geom_blank(data=axis_limits,aes(x=Group,y=Limit),inherit.aes=FALSE)+
  scale_y_continuous(breaks=seq(-8,8,by=2),expand=expansion(mult=.025))+
  geom_hline(yintercept=0,linetype='dashed',linewidth=.35,colour='grey50')+
  geom_boxplot(aes(group=Group),outlier.shape=NA,fill='grey95',width=.6,linewidth=.45)+
  geom_point(aes(colour=Group,shape=Group),size=1.5,alpha=.8,
   position=position_jitter(width=.1,height=0,seed=20261109))+
  geom_text(data=undefined,aes(x=Group,y=0,label='ES not estimable'),inherit.aes=FALSE,
   colour='grey35',size=3.5)+
  facet_wrap(~Trait_label,ncol=2,scales='free_y',drop=FALSE)+
  scale_colour_manual(name='Regime / pressure',values=colours,labels=labels,drop=FALSE)+
  scale_shape_manual(name='Regime / pressure',values=shapes,labels=labels,drop=FALSE)+
  scale_x_discrete(drop=FALSE,labels=c('P Low\nn = 30','P High\nn = 16','I Low\nn = 14','I High\nn = 13'))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
  labs(title=unname(titles[family]),subtitle='CF v3 | 73 unaveraged libraries | fixed verified publication null',
   x=NULL,y='Standardized effect size (ES)')+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
   strip.text=element_text(size=10),legend.position='bottom',legend.text=element_text(size=10),
   legend.title=element_text(size=10),axis.text.x=element_text(size=8),plot.margin=margin(10,12,10,12))
 height <- max(4.8,1.7*ceiling(length(traits)/2)+1.6)
 filename <- sprintf('ES_S%d_%s_COLOURED_FIXED_SCALE',i,family)
 ggsave(file.path(out,paste0(filename,'.png')),p,width=9,height=height,dpi=400,bg='white')
 ggsave(file.path(out,paste0(filename,'.tiff')),p,width=9,height=height,dpi=600,compression='lzw',bg='white')
 ggsave(file.path(out,paste0(filename,'.pdf')),p,width=9,height=height,bg='white')
 captions <- c(captions,paste0('ES Figure S',i,'. ',unname(titles[family]),
 '. Trait-specific standardized effect sizes for all 73 unaveraged sequencing libraries using CF v3-corrected relative abundance and the verified 9,999-permutation null. Colours and symbols identify the four regime-pressure groups. P denotes perennial and I intermittent. Boxes show medians and interquartile ranges; whiskers extend to 1.5 IQR. Points are libraries, including technical replicates. Y-axes are symmetric about zero, with ticks every two ES units and a minimum range of -3 to 3; limits are extended only where needed to retain every point. The dashed line marks the null expectation, not a test of group means against zero. ',
 if(family=='Colony')'Stellate_colony and Other_colony have zero null variance; their panels are labelled not estimable and contain no artificial zero ES values.'else'',
 ' Figure identifiers are local to this ES collection and do not replace existing supplement numbering.'))
}
writeLines(captions,file.path(out,'ES_trait_group_captions.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))

# Figure 5 retains its coloured groups and clustered confidence intervals.
library(sandwich)
z <- d[d$Trait=='Not_attached'&is.finite(d$ES),]
stopifnot(nrow(z)==73)
z$Field_station <- meta$field_station[match(z$Unit_ID,meta$sample_id)]
fit <- lm(ES~0+Group,z)
V <- vcovCL(fit,cluster=z$Field_station,type='HC1')
em <- data.frame(Group=levels(z$Group),Mean=coef(fit),SE=sqrt(diag(V)))
em$Lower95 <- em$Mean-qt(.975,37)*em$SE
em$Upper95 <- em$Mean+qt(.975,37)*em$SE
p5 <- ggplot(z,aes(Group,ES))+geom_hline(yintercept=0,linetype='dashed',colour='grey50')+
 geom_point(aes(colour=Group,shape=Group),size=2.7,alpha=.8,
 position=position_jitter(width=.1,height=0,seed=20261007))+
 geom_errorbar(data=em,aes(Group,ymin=Lower95,ymax=Upper95),inherit.aes=FALSE,width=.14,linewidth=.7)+
 geom_point(data=em,aes(Group,Mean),inherit.aes=FALSE,shape=23,fill='white',size=3.7)+
 scale_colour_manual(name='Regime / pressure',values=colours,labels=labels)+
 scale_shape_manual(name='Regime / pressure',values=shapes,labels=labels)+
 scale_x_discrete(labels=c('Perennial\nLow','Perennial\nHigh','Intermittent\nLow','Intermittent\nHigh'))+
 scale_y_continuous(limits=c(-3,3),breaks=c(-2,0,2),expand=expansion(mult=.025))+
 guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
 labs(x=NULL,y='Standardized effect size (ES)',title='Not attached',
 subtitle='CF v3 | 73 libraries | pressure q: OLS 0.0010; cluster 0.0417')+
 theme_classic(base_size=13)+theme(legend.position='bottom',legend.text=element_text(size=11),
 legend.title=element_text(size=11),plot.margin=margin(12,14,12,14))
ggsave(file.path(out,'Figure_5_ES_COLOURED_FIXED_SCALE.png'),p5,width=8.8,height=6.8,dpi=600,bg='white')
ggsave(file.path(out,'Figure_5_ES_COLOURED_FIXED_SCALE.tiff'),p5,width=8.8,height=6.8,dpi=600,compression='lzw',bg='white')
ggsave(file.path(out,'Figure_5_ES_COLOURED_FIXED_SCALE.pdf'),p5,width=8.8,height=6.8,bg='white')
