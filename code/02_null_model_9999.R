# FULL NULL MODEL FROM THE BASIC TABLES -- no custom/anonymous functions.
# Run from the repository root: source('code/02_null_model_9999.R')
# Native R RNG differs from the original NumPy PCG64 stream; same null method,
# not the identical historical random draws. Published reference files stay intact.
stopifnot(file.exists('data/taxa_trait.csv'))
nr <- 9999L
seed <- 20261007L
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(seed)
out_null <- 'null_rerun_R'
dir.create(out_null,showWarnings=FALSE)
trait <- read.csv('data/taxa_trait.csv',check.names=FALSE,encoding='UTF-8')
raw <- read.csv('data/sample_taxa_rawread.csv',check.names=FALSE,encoding='UTF-8')
env <- read.csv('data/sample_env.csv',check.names=FALSE,encoding='UTF-8')
dictionary <- read.csv('data/trait_dictionary.csv',check.names=FALSE,encoding='UTF-8')
stopifnot(nrow(raw)==73,nrow(env)==73,!anyDuplicated(raw$sample_id),
          !anyDuplicated(trait$taxon),identical(raw$sample_id,env$sample_id))
taxa_all <- names(raw)[-1]
idx <- match(taxa_all,trait$taxon)
stopifnot(!anyNA(idx))
read_matrix <- as.matrix(raw[,-1])
rownames(read_matrix) <- raw$sample_id
stopifnot(all(is.finite(read_matrix)),all(read_matrix>=0),all(rowSums(read_matrix)>0))
factor_cf <- trait$CF_v3[idx]
covered <- is.finite(factor_cf)&factor_cf>0
# Derive corrected reads directly from raw reads and the authoritative CF field.
corrected <- sweep(read_matrix[,covered,drop=FALSE],2,factor_cf[covered],'/')
record <- list()
cover_record <- list()
started <- Sys.time()
for(weight in c('Raw_all','Raw_matched','CF3_matched')) {
 if(weight=='Raw_all') {
  abundance <- read_matrix;taxa <- taxa_all
 } else if(weight=='Raw_matched') {
  abundance <- read_matrix[,covered,drop=FALSE];taxa <- taxa_all[covered]
 } else {
  abundance <- corrected;taxa <- taxa_all[covered]
 }
 # Retain zero-read taxa in the annotation pool, as in the original final null.
 relative <- abundance/rowSums(abundance)
 q <- as.matrix(trait[match(taxa,trait$taxon),dictionary$Column])
 stopifnot(all(is.na(q)|q==0|q==1))
 for(tr in dictionary$Column) {
  known <- !is.na(q[,tr])
  labels <- q[known,tr]
  known_weight <- rowSums(relative[,known,drop=FALSE])
  stopifnot(length(labels)>0,all(known_weight>0))
  normalized <- relative[,known,drop=FALSE]/known_weight
  observed_p <- as.numeric(normalized%*%labels)
  observed <- 1/(observed_p^2+(1-observed_p)^2)
  # One independent permutation per trait/weighting, shared by all 73 samples.
  # Explicit loop is intentional: there is no function definition or apply closure.
  permuted <- matrix(NA_real_,nrow=length(labels),ncol=nr)
  for(permutation in seq_len(nr)) {
   permuted[,permutation] <- labels[sample.int(length(labels),replace=FALSE)]
  }
  p_null <- normalized%*%permuted
  diversity_null <- 1/(p_null^2+(1-p_null)^2)
  mean_null <- rowMeans(diversity_null)
  centred <- diversity_null-mean_null
  sd_null <- sqrt(rowSums(centred^2)/(nr-1))
  # Remove round-off variance for mathematically constant null distributions.
  spread_null <- apply(diversity_null,1,max)-apply(diversity_null,1,min)
  sd_null[spread_null<1e-13] <- 0
  if(length(unique(labels))==1) {
   mean_null[] <- 1;sd_null[] <- 0
  }
  record[[length(record)+1L]] <- data.frame(Unit_ID=raw$sample_id,Level='Library73',
   Weighting=weight,Trait=tr,Observed=observed,Null_mean=mean_null,Null_sd=sd_null,N_random=nr)
  cover_record[[length(cover_record)+1L]] <- data.frame(Unit_ID=raw$sample_id,
   Weighting=weight,Trait=tr,Known_taxa=sum(known),Known_weight=known_weight)
  cat(weight,tr,'completed;',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'s\n')
  flush.console()
 }
}
null_rerun <- do.call(rbind,record)
coverage_rerun <- do.call(rbind,cover_record)
stopifnot(nrow(null_rerun)==73*42*3,all(is.finite(null_rerun$Observed)),
 all(is.finite(null_rerun$Null_mean)),all(is.finite(null_rerun$Null_sd)))
write.csv(null_rerun,file.path(out_null,'null_summary_9999.csv'),row.names=FALSE)
write.csv(coverage_rerun,file.path(out_null,'known_trait_coverage.csv'),row.names=FALSE)
writeLines(c(paste('Seed:',seed),paste('Permutations:',nr),paste('RNG:',paste(RNGkind(),collapse=', ')),
 'Permutation unit: known taxon labels; same permutation shared by all 73 samples.',
 'Pool: all annotated taxa per weighting, including zero-read taxa.',
 'ES: (observed diversity - mean null diversity) / sample SD null diversity.',
 'Original publication used NumPy PCG64. Native R random draws differ.',
 paste('Started:',started),paste('Finished:',Sys.time()),capture.output(sessionInfo())),
 file.path(out_null,'provenance.txt'))
# Continue with the same final models/figures, now using the freshly simulated null.
options(cyprus.null_source='R_rerun')
source('code/01_ES_final.R')
