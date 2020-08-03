###GI frequency
ps <- readRDS("data/ps_noDuplicates.RDS")
sample_data(ps)$bowel_norm <- ifelse(sample_data(ps)$bowel.freq == 1, "Normal", "Abnormal")
ps_aut <- subset_samples(ps, Treatment == "Aut")
ps_nt <- subset_samples(ps, Treatment == "Control")


#Frequency
par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
pie(table(sample_data(ps_aut)$bowel.freq) / nsamples(ps_aut), labels = c("0", "1", "1.5", "2", "3", "4"),
    sub = "ASD")
pie(table(sample_data(ps_nt)$bowel.freq) / nsamples(ps_nt), labels = c("0", "1", "1.5", "2", "3", "4"),
    sub = "NT")
mtext("Typical number of bowel movements a day" ,outer = T, cex = 1.5)


#Quality
par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
pie(table(sample_data(ps_aut)$bowel.qual) / nsamples(ps_aut), labels = c("Constipation", "Diarrhea", "Normal", "Unknown"),
    sub = "ASD")
pie(table(sample_data(ps_nt)$bowel.qual) / nsamples(ps_nt), labels = c("Constipation", "Diarrhea", "Normal", "Unknown"),
    sub = "NT")
mtext("Typical quality of bowel movements" ,outer = T, cex = 1.5)


#Subdivide into low and high frequency bowel movement
ps_aut_lowFreq <- subset_samples(ps_aut, bowel.freq %in% c(0,1,1.5))
ps_aut_highFreq <- subset_samples(ps_aut, bowel.freq %in% c(2,3,4))

ps_nt_lowFreq <- subset_samples(ps_nt, bowel.freq %in% c(0,1,1.5))
ps_nt_highFreq <- subset_samples(ps_nt, bowel.freq %in% c(2,3,4))


####################################################################
#TODO: clean script so you can import properly
#source('ASD_microbiome_final.R')
runDESeq <- function(ps, confounders = NA){
  if(is.na(confounders)){
    form <- as.formula('~Treatment')
  }else{
    form <- as.formula(paste('~Treatment', confounders, sep = "+"))
  }
  print(form)
  diagdds = phyloseq_to_deseq2(ps, form )
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("Treatment", "Aut", "Control"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

run_metagenom_seq<-function(ps,maxit, confounders){
  p_metag<-phyloseq_to_metagenomeSeq(ps)
  #filtering at least 4 samples 
  p_metag= cumNorm(p_metag, p=0.75)
  normFactor =normFactors(p_metag)
  normFactor =log2(normFactor/median(normFactor) + 1)
  #mod = model.matrix(~ASDorNeuroT +PairASD+ normFactor)
  
  if(is.na(confounders)){
    form <- as.formula('~Treatment + normFactor')
  }else{
    form <- as.formula(paste('~Treatment + normFactor', confounders, sep = "+"))
  }
  print(form)
  mod = model.matrix(form, data = pData(p_metag))
  settings =zigControl(maxit =maxit, verbose = TRUE)
  #settings =zigControl(tol = 1e-5, maxit = 30, verbose = TRUE, pvalMethod = 'bootstrap')
  fit =fitZig(obj = p_metag, mod = mod, useCSSoffset = FALSE, control = settings)
  #then finding which tax have to have #38 -> 116: mean 51
  res_fit<-MRtable(fit, number = length(fit$taxa))
  res_fit <- res_fit[!is.na(res_fit$adjPvalues), ]
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

ps <- readRDS("data/ps_noDuplicates.RDS")

sample_data(ps)$bowelFreq_Cat <- ifelse(sample_data(ps)$bowel.freq %in% c(0,1,1.5), "low", "high")
levels(sample_data(ps)$bowel.freq) <- c(-1, levels(sample_data(ps)$bowel.freq))
sample_data(ps)$bowel.freq[is.na(sample_data(ps)$bowel.freq)] <- -1
sample_data(ps)$bowel.freq <- as.numeric(sample_data(ps)$bowel.freq)

ps <- removeASD_noDiag(ps)
###Filter and normalize
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
ps_filt_norm <- deSeqNorm(ps_filt)



###Run DESEQ
deseq_res <- runDESeq(ps_filt, confounders = "bowel.freq")
deseq_res$enriched <- ifelse(deseq_res$log2FoldChange > 0,"Aut", "Control")


###Run ZIG model fitting and prediction
zig_res <- run_metagenom_seq(ps_filt,30, confounders = c("bowel.freq"))
zig_res <- data.frame(cbind(zig_res), tax_table(ps_filt)[rownames(zig_res), ])
zig_res$enriched <- ifelse(zig_res$TreatmentControl < 0, "Aut", "Control")



