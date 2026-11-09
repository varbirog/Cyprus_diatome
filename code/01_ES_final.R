# Final ES package, all 73 libraries; no station averaging or manuscript edits.
# Uses the verified 9,999-permutation null summaries; does not rerun the null.
# Set the working directory to this repository root before running.
stopifnot(file.exists('data/sample_env.csv'))
library(ggplot2)
library(sandwich)
library(car)
b <- 'data'
null_source <- getOption('cyprus.null_source','published')
stopifnot(null_source %in% c('published','R_rerun'))
out <- if(null_source=='published') 'results' else 'results_R_null'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'figures'),showWarnings=FALSE)
meta <- read.csv(file.path(b,'sample_env.csv'),check.names=FALSE,encoding='UTF-8')
meta$Sample_ID <- meta$sample_id
meta$Station <- meta$field_station
null_file <- if(null_source=='published') file.path(b,'null_summary_9999.csv') else 'null_rerun_R/null_summary_9999.csv'
null <- read.csv(null_file)
null <- null[null$Level=='Library73',]
stopifnot(nrow(meta)==73,nrow(null)==73*42*3,all(null$N_random==9999))
null$ES <- (null$Observed-null$Null_mean)/null$Null_sd
null$ES[!is.finite(null$ES)|null$Null_sd<=0] <- NA_real_
# Verify all observed diversity values directly from raw reads and final traits.
trait <- read.csv('data/taxa_trait.csv',check.names=FALSE,encoding='UTF-8')
raw <- read.csv('data/sample_taxa_rawread.csv',check.names=FALSE,encoding='UTF-8')
cf3 <- read.csv('data/sample_taxa_cf3read.csv',check.names=FALSE,encoding='UTF-8')
dictionary <- read.csv('data/trait_dictionary.csv',check.names=FALSE,encoding='UTF-8')
stopifnot(identical(raw$sample_id,meta$Sample_ID),identical(cf3$sample_id,meta$Sample_ID))
for(w in unique(null$Weighting)) {
 if(w=='Raw_all') {
  abund <- as.matrix(raw[,-1]);taxa <- names(raw)[-1]
 } else if(w=='Raw_matched') {
  taxa <- names(cf3)[-1];abund <- as.matrix(raw[,taxa])
 } else {
  taxa <- names(cf3)[-1];abund <- as.matrix(cf3[,-1])
 }
 q <- as.matrix(trait[match(taxa,trait$taxon),dictionary$Column])
 stopifnot(!anyNA(match(taxa,trait$taxon)))
 rel <- abund/rowSums(abund)
 known <- !is.na(q)
 q[is.na(q)] <- 0
 prop <- (rel%*%q)/(rel%*%(known*1))
 observed <- 1/(prop^2+(1-prop)^2)
 for(tr in dictionary$Column) {
  ix <- null$Weighting==w&null$Trait==tr
  stopifnot(max(abs(null$Observed[ix]-observed[match(null$Unit_ID[ix],meta$Sample_ID),tr]))<1e-10)
 }
}
cat('All 9198 observed diversity values independently verified from the basic tables.\n')

previous <- read.csv('reference/ES_verified.csv')
key <- paste(null$Unit_ID,null$Weighting,null$Trait)
oldkey <- paste(previous$Unit_ID,previous$Weighting,previous$Trait)
if(null_source=='published') {
 stopifnot(max(abs(null$ES-previous$ES[match(key,oldkey)]),na.rm=TRUE)<1e-10)
}
null$Field_station <- meta$Station[match(null$Unit_ID,meta$Sample_ID)]
null$regime <- factor(meta$regime[match(null$Unit_ID,meta$Sample_ID)],levels=c('Perennial','Intermittent'))
null$pressure <- factor(meta$pressure[match(null$Unit_ID,meta$Sample_ID)],levels=c('Low','High'))
rows <- list()
for(w in unique(null$Weighting))for(tr in unique(null$Trait)) {
 z <- null[null$Weighting==w&null$Trait==tr&is.finite(null$ES),]
 if(nrow(z)==0)next
 stopifnot(nrow(z)==73)
 fit <- lm(ES~regime*pressure,z,contrasts=list(regime=contr.sum,pressure=contr.sum))
 ols <- as.data.frame(Anova(fit,type=3))
 se_hc3 <- sqrt(diag(vcovHC(fit,type='HC3')))
 se_cluster <- sqrt(diag(vcovCL(fit,cluster=z$Field_station,type='HC1')))
 for(j in 2:4) {
  term <- c('regime','pressure','regime:pressure')[j-1]
  rows[[length(rows)+1]] <- data.frame(Weighting=w,Trait=tr,Term=term,N=73,
    Coefficient=coef(fit)[j],OLS_F=ols[term,'F value'],OLS_P=ols[term,'Pr(>F)'],
    HC3_SE=se_hc3[j],HC3_P=2*pt(-abs(coef(fit)[j]/se_hc3[j]),69),
    Cluster_SE=se_cluster[j],Cluster_P=2*pt(-abs(coef(fit)[j]/se_cluster[j]),37))
 }
}
models <- do.call(rbind,rows)
models$OLS_Q <- models$HC3_Q <- models$Cluster_Q <- NA_real_
for(w in unique(models$Weighting))for(term in unique(models$Term)) {
 ix <- models$Weighting==w&models$Term==term
 models$OLS_Q[ix] <- p.adjust(models$OLS_P[ix],'BH',n=42)
 models$HC3_Q[ix] <- p.adjust(models$HC3_P[ix],'BH',n=42)
 models$Cluster_Q[ix] <- p.adjust(models$Cluster_P[ix],'BH',n=42)
}
old <- read.csv('reference/ES_cluster_verified.csv')
old <- old[is.finite(old$FDR),]
mkey <- paste(models$Weighting,models$Trait,models$Term)
if(null_source=='published') {
 stopifnot(max(abs(models$Cluster_Q-old$FDR[match(mkey,paste(old$Weighting,old$Trait,old$Term))]))<1e-10)
}
old_ols <- read.csv('reference/ES_OLS_verified.csv')
if(null_source=='published') {
 stopifnot(max(abs(models$OLS_Q-old_ols$FDR[match(mkey,paste(old_ols$Weighting,old_ols$Trait,old_ols$Term))]))<1e-10)
}
write.csv(null,file.path(out,'S7_ES_all_73_profiles.csv'),row.names=FALSE,na='')
write.csv(models,file.path(out,'S8_ES_models_OLS_HC3_cluster.csv'),row.names=FALSE,na='')
write.csv(models[models$OLS_Q<.05,],file.path(out,'S7_ES_significant_OLS.csv'),row.names=FALSE)
write.csv(models[models$Cluster_Q<.05,],file.path(out,'S8_ES_significant_cluster.csv'),row.names=FALSE)
means <- aggregate(ES~Weighting+Trait+regime+pressure,null,mean)
write.csv(means,file.path(out,'S7_ES_group_means.csv'),row.names=FALSE)
# Figure 4: continuous ES colours indicate sign, not categorical group colours.
h <- means[means$Weighting=='CF3_matched',]
dictionary <- read.csv(file.path(b,'trait_dictionary.csv'))
variable_traits <- dictionary$Column[dictionary$Column %in% h$Trait]
h$Trait <- factor(gsub('_',' ',h$Trait),levels=rev(gsub('_',' ',variable_traits)))
h$Group <- factor(paste(h$regime,h$pressure,sep=' / '),
  levels=c('Perennial / Low','Perennial / High','Intermittent / Low','Intermittent / High'))
p4 <- ggplot(h,aes(Group,Trait,fill=ES))+geom_tile(colour='white',linewidth=.25)+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,name='Mean ES')+
 labs(x=NULL,y=NULL,subtitle='CF v3 | 73 unaveraged libraries | 40 estimable trait categories')+
 theme_minimal(base_size=11)+theme(panel.grid=element_blank(),axis.text.y=element_text(size=9),
 axis.text.x=element_text(angle=25,hjust=1),plot.margin=margin(12,12,12,12))
# Figure 5: preserve all observations and use station-clustered confidence intervals.
z <- null[null$Weighting=='CF3_matched'&null$Trait=='Not_attached',]
z$Group <- factor(paste(z$regime,z$pressure,sep='_'),
 levels=c('Perennial_Low','Perennial_High','Intermittent_Low','Intermittent_High'))
fit <- lm(ES~0+Group,z)
V <- vcovCL(fit,cluster=z$Field_station,type='HC1')
em <- data.frame(Group=levels(z$Group),Mean=coef(fit),SE=sqrt(diag(V)))
em$Lower95 <- em$Mean-qt(.975,37)*em$SE
em$Upper95 <- em$Mean+qt(.975,37)*em$SE
write.csv(em,file.path(out,'S7_Not_attached_clustered_model_means.csv'),row.names=FALSE)
colours <- c(Perennial_Low='#0072B2',Perennial_High='#56B4E9',Intermittent_Low='#D55E00',Intermittent_High='#E69F00')
shapes <- c(Perennial_Low=16,Perennial_High=17,Intermittent_Low=15,Intermittent_High=3)
labels <- c('Perennial / Low','Perennial / High','Intermittent / Low','Intermittent / High')
m <- models[models$Weighting=='CF3_matched'&models$Trait=='Not_attached'&models$Term=='pressure',]
p5 <- ggplot(z,aes(Group,ES))+geom_hline(yintercept=0,linetype='dashed',colour='grey50')+
 geom_point(aes(colour=Group,shape=Group),size=2.7,alpha=.8,
 position=position_jitter(width=.1,height=0,seed=20261007))+
 geom_errorbar(data=em,aes(Group,y=Mean,ymin=Lower95,ymax=Upper95),inherit.aes=FALSE,width=.14,linewidth=.7)+
 geom_point(data=em,aes(Group,Mean),inherit.aes=FALSE,shape=23,fill='white',size=3.7)+
 scale_colour_manual(name='Regime / pressure',values=colours,labels=labels)+
 scale_shape_manual(name='Regime / pressure',values=shapes,labels=labels)+
 scale_x_discrete(labels=c('Perennial\nLow','Perennial\nHigh','Intermittent\nLow','Intermittent\nHigh'))+
 guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
 labs(x=NULL,y='Standardized effect size (ES)',title='Not attached',
 subtitle=sprintf('CF v3 | 73 libraries | pressure q: OLS %.4f; cluster %.4f',m$OLS_Q,m$Cluster_Q))+
 theme_classic(base_size=13)+theme(legend.position='bottom',legend.text=element_text(size=11),
 legend.title=element_text(size=11),plot.margin=margin(12,14,12,14))
ggsave(file.path(out,'figures/Figure_4_ES_FINAL.png'),p4,width=8.4,height=10,dpi=600,bg='white')
ggsave(file.path(out,'figures/Figure_4_ES_FINAL.tiff'),p4,width=8.4,height=10,dpi=600,compression='lzw',bg='white')
ggsave(file.path(out,'figures/Figure_4_ES_FINAL.pdf'),p4,width=8.4,height=10,bg='white')
ggsave(file.path(out,'figures/Figure_5_ES_FINAL.png'),p5,width=8.8,height=6.8,dpi=600,bg='white')
ggsave(file.path(out,'figures/Figure_5_ES_FINAL.tiff'),p5,width=8.8,height=6.8,dpi=600,compression='lzw',bg='white')
ggsave(file.path(out,'figures/Figure_5_ES_FINAL.pdf'),p5,width=8.8,height=6.8,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Completed ES analysis with null source:',null_source,'\n')
print(models[models$Weighting=='CF3_matched'&(models$OLS_Q<.05|models$Cluster_Q<.05),])
