library(metagenomeSeq)
library(phyloseq)

#setwd("/Users/mmdavid/Documents/microbiome_paper_FINAL/ASD_microbiome/")
#ok working with metagenomeseq now
set.seed(100)
source ("scripts/getTables.R")


#setting up the sample to use
find_zero<-function(ZIGObs,phyloOb){
  tmp_Ob<-otu_table(prune_taxa(rownames(ZIGObs),phyloOb))
  apply(1,function(x) {length(x[x!=0])})
}


filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold = percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(otu_table(ps), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  otu_table(ps) <- otu_table(ps)[toKeep,]
  return(ps)
}

run_metagenom_seq<-function(ps,maxit){
  p_metag<-phyloseq_to_metagenomeSeq(ps)
  #filtering at least 4 samples 
  p_metag= cumNorm(p_metag, p=0.75)
  normFactor =normFactors(p_metag)
  normFactor =log2(normFactor/median(normFactor) + 1)
  #mod = model.matrix(~ASDorNeuroT +PairASD+ normFactor)
  mod = model.matrix(~Treatment + normFactor, data = pData(p_metag))
  settings =zigControl(maxit =maxit, verbose = TRUE)
  #settings =zigControl(tol = 1e-5, maxit = 30, verbose = TRUE, pvalMethod = 'bootstrap')
  fit =fitZig(obj = p_metag, mod = mod, useCSSoffset = FALSE, control = settings)
  #then finding which tax have to have #38 -> 116: mean 51
  res_fit<-MRtable(fit, number = length(fit$taxa))
  res_fit<-res_fit[res_fit$adjPvalues<0.05,]
  #finally remove the ones that are not with enough samples
  #mean_sample<-mean(calculateEffectiveSamples(fit))
  #res_fit<-res_fit[res_fit$`counts in group 0` & res_fit$`counts in group 1` > mean_sample,]
  Min_effec_samp<-calculateEffectiveSamples(fit)
  Min_effec_samp<-Min_effec_samp[ names(Min_effec_samp)  %in% rownames(res_fit)]
  res_fit$Min_sample<-Min_effec_samp
  res_fit<-res_fit[res_fit$`+samples in group 0` >= Min_effec_samp & res_fit$`+samples in group 1` >= Min_effec_samp,]
  return(res_fit)
}


main <- function(){
  #setwd("~/Documents/Microbiome_ASD_16S")
  ps <- readRDS("data/ps_filt_norm_117.RDS")
  ps_03 <- filterTaxaByPrevolence(ps, .03)
  res <- run_metagenom_seq(ps_03,30 )
  #Group 0 is aut
  
  taxa <- as.data.frame(tax_table(ps)[rownames(res), ])
  taxa$adjPval <- res$adjPvalues
  taxa$enriched <- ifelse(res$TreatmentControl < 0, "Aut", "Control")
  
  saveRDS(taxa, "metagenomseq/filter03_117s.rds")
}
 main()
