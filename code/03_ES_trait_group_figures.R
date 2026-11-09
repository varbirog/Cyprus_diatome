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
out <- 'ES_trait_group_figures'
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
 p <- ggplot(usable,aes(Group,ES))+
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
 filename <- sprintf('ES_S%d_%s_FINAL',i,family)
 ggsave(file.path(out,paste0(filename,'.png')),p,width=9,height=height,dpi=400,bg='white')
 ggsave(file.path(out,paste0(filename,'.tiff')),p,width=9,height=height,dpi=600,compression='lzw',bg='white')
 ggsave(file.path(out,paste0(filename,'.pdf')),p,width=9,height=height,bg='white')
 captions <- c(captions,paste0('ES Figure S',i,'. ',unname(titles[family]),
 '. Trait-specific standardized effect sizes for all 73 unaveraged sequencing libraries using CF v3-corrected relative abundance and the verified 9,999-permutation null. Colours and symbols identify the four regime-pressure groups. P denotes perennial and I intermittent. Boxes show medians and interquartile ranges; whiskers extend to 1.5 IQR. Points are libraries, including technical replicates. Vertical scales vary among traits. The dashed line marks the null expectation, not a test of group means against zero. ',
 if(family=='Colony')'Stellate_colony and Other_colony have zero null variance; their panels are labelled not estimable and contain no artificial zero ES values.'else'',
 ' Figure identifiers are local to this ES collection and do not replace existing supplement numbering.'))
}
writeLines(captions,file.path(out,'ES_trait_group_captions.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
