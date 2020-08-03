#dsetwd("C:/Users/Christine/Documents/Lab/Microbiome_ASD_16S/")
loc.mapping = "mapping/mapping_final_withOutliers.txt" ##This mapping file does not have the smaples with low quality


##################################################################################################
########################## Pairwise analysis using OTU
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold = percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(otu_table(ps), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  otu_table(ps) <- otu_table(ps)[toKeep,]
  return(ps)
}

deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ Treatment )
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund))
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table(ps))
  return(ps_deSeq)
}

getDifferenceMeans <- function(table, mapping){
  pairs_otu <- splitAutControl(table,mapping)
  diffPairs_otus<-pairs_otu[[1]]-pairs_otu[[2]]
  diffPairs_otus<-diffPairs_otus[rowSums(diffPairs_otus) != 0,] 
  return(apply(diffPairs_otus,1,mean))
}

permutation <- function(table, mapping){
  pairs_otu <- splitAutControl3_shuff(table,mapping)
  #tmp<-min(length(pairs_otu[[1]]),length(pairs_otu[[2]]))
  diffPairs_otus<-pairs_otu[[1]]-pairs_otu[[2]]
  return(apply(diffPairs_otus,1,mean))
  #return (diffPairs_otus)
}

runSimulations <- function(table, mapping, numSims){
  null_difference_means <- matrix(ncol = nrow(table), nrow=numSims)
  colnames(null_difference_means) <- rownames(table)
  for (i in 1:numSims){
    cat (i,"\t")
    null_difference_means[i,] <- as.numeric(permutation(table, mapping))
  }
  return(null_difference_means)
}

calcPVals <- function(null_difference_means, actual_difference_means){
  #ok now the test: using the area under the curve 
  p.vals = c()
  enrichments <- c()
  for (i in 1: ncol(null_difference_means)){ # i is an otu index
    relevant = c()
    enriched <- c()
    null_dist <- null_difference_means[ , i]
    actual_value <- actual_difference_means[[i]]
    if(actual_value > mean(null_dist)){
      pval = (1+sum(null_dist > actual_value)) / 10001
      enriched <- "Aut"}
    else{
      pval = (1 + sum(null_dist < actual_value)) / 10001
      enriched <- "Control"
    }
    pval = pval * 2
    p.vals <-c(p.vals , pval)
    enrichments <- c(enrichments, enriched)
  }
  return(list(p.vals, enrichments))
}

##############################################
##############################################
##############################################
##############################################
##############################################NOT WORKING WITH PRELOADED RESULTS
mcPermutationTest <- function(ps, null_difference_means = NA, numSims = 10000){
  #load data
  #load("data/phyloseq_obj_deseq")
  ps_use = ps

  mapping_otus <- as.data.frame(sample_data(ps_use))
  otu_table <- otu_table(ps_use)
  
  #OK let's do the pair analysis, here we're separating by pairs: total 62 pairs, because we have fmailies with 3 kids
  
  if(all(is.na(null_difference_means))){
    null_difference_means <- runSimulations(otu_table, mapping_otus, numSims) #returns a simulation X factor table
  }
  
  
  actual_difference_means <- getDifferenceMeans(otu_table, mapping_otus)
  
  #the function "getDifferenceMeans" removes the factors that have zeros only: 
  #this also ensures the columns are in the same order
  null_difference_means <- null_difference_means[,names(actual_difference_means)]
  colnames(null_difference_means) == names(actual_difference_means)
  #saveRDS(null_difference_means, paste("pair_analysis_OTU_picrust/nullDistributionShapes/null_dist_10000_", dataType, "_135.rds", sep = ""))
  
  tmp <- calcPVals(null_difference_means, actual_difference_means)
  pvals <- tmp[[1]]
  enrichments <- tmp[[2]]
  pvals_adjusted <- p.adjust(pvals, method = "fdr")
  
  #Adjust with tree
  #names(pvals) <- colnames(null_difference_means)
  #edge_list <- get.edgelist(as.igraph(phy_tree(ps_use)))
  #hfdr_res <- hFDR.adjust(pvals, edge_list, .05)
  #summary(hfdr_res)
  
  df <- data.frame(pvals = pvals, pvals_adj = pvals_adjusted, enriched = enrichments,
                   seqs = colnames(null_difference_means))
  df_sig <- df[df$pvals_adj < .1, ]
  tax = as.data.frame(tax_table(ps_use))[as.character(df_sig$seqs), ]
  return(cbind(df_sig, tax))
}

##############################################
##############################################
##############################################
##############################################
##############################################

#Use to run independently
main <- function(){
  ps <- readRDS("data/ps_filt_norm_117.RDS")
  ps <- subset_samples(ps, !SampleID %in% c("45.1", "233")) #take out the third siblings
  pairs_to_keep <- sample_data(ps)$Pair[duplicated(sample_data(ps)$Pair)]
  ps <- subset_samples(ps, Pair %in% pairs_to_keep)
  
  ps <- filterTaxaByPrevolence(ps, .03)
  ps_deSeq <- deSeqNorm(ps)
  
  
  #ps_use <- ps_deSeq
  res = mcPermutationTest(ps_deSeq, numSims = 10000)
  res$nb_non_zero = apply(otu_table(ps)[res$seqs, ], 1, function(tax) return(sum(tax>0)))
  outfile = "pair_analysis_OTU_picrust/mc_03filter_110s_deseq_adj.RDS"
  saveRDS(res, file=outfile)
}


##############################################
##############################################
##############################################
##############################################
##############################################
#Results using this command below
run <- function(ps, norm = "", numSims = 100){
  #load data
  #load("data/phyloseq_obj_deseq")
  ps_use = ps
  
  mapping_otus <- as.data.frame(sample_data(ps_use))
  otu_table <- otu_table(ps_use)
  
  #OK let's do the pair analysis, here we're separating by pairs: total 62 pairs, because we have fmailies with 3 kids
  null_difference_means <- runSimulations(otu_table, mapping_otus, numSims) #returns a simulation X factor table
  actual_difference_means <- getDifferenceMeans(otu_table, mapping_otus)
  
  #the function "getDifferenceMeans" removes the factors that have zeros only: 
  #this also ensures the columns are in the same order
  null_difference_means <- null_difference_means[,names(actual_difference_means)]
  colnames(null_difference_means) == names(actual_difference_means)
  #saveRDS(null_difference_means, paste("pair_analysis_OTU_picrust/nullDistributionShapes/null_dist_10000_", dataType, "_135.rds", sep = ""))
  
  tmp <- calcPVals(null_difference_means, actual_difference_means)
  pvals <- tmp[[1]]
  enrichments <- tmp[[2]]
  pvals_adjusted <- p.adjust(pvals, method = "fdr")
  
  df <- data.frame(pvals = pvals, pvals_adj = pvals_adjusted, enriched = enrichments,
                   seqs = colnames(null_difference_means))
  df_sig <- df[df$pvals_adj < .1, ]
  tax = as.data.frame(tax_table(ps_use))[as.character(df_sig$seqs), ]
  return(cbind(df_sig, tax))
}







##################################################
numSimsNecessary <- function(){
  pvalues <- list()
  MAPPING <- sample_data(getMapping(loc.mapping))
  MAPPING <- MAPPING[!MAPPING$SampleID %in% c('367', '173','168.1' ,'384', '163', '211'),] # 173 and 367 are extreme outliers after CSSnorm, 
  
  OTU_ALL <- otu_table(phyloseq_obj_all)
  OTU_ALL <- OTU_ALL[,colnames(OTU_ALL) %in% MAPPING$SampleID]
  totalSims <- 50000
  batchSize <- 500
  
  for(iteration in seq(1, totalSims / batchSize)){
    null_difference_means_batch <- runSimulations(OTU_ALL, MAPPING, batchSize) #each column is an otu index, rows are simulations
    if(iteration == 1){
      null_difference_means_current <- null_difference_means_batch
      null_difference_means_prev <- matrix(rep(0, nrow(null_difference_means_current) * ncol(null_difference_means_current)), nrow = nrow(null_difference_means_current))
    }else{
      null_difference_means_current <- rbind(null_difference_means_prev, null_difference_means_batch)
    }
    asList_current <- lapply(seq(1:ncol(null_difference_means_current)), function(i) null_difference_means_current[,i])
    asList_prev <- lapply(seq_len(ncol(null_difference_means_prev)), function(i) null_difference_means_prev[,i])
    
    taxa_delta_ks <- mapply(ks.test, asList_current, asList_prev)
    # taxa_delta_ks_normal <- lapply(asList_current, function(taxa_dist) return(shapiro.test(taxa_dist)$p.value))
    
    pvalues[[iteration]] <- unlist(taxa_delta_ks[seq(2, length(taxa_delta_ks), by = 5)])
    # pvalues_normal[[iterations]] <-unlist(taxa_delta_ks_normal)
    null_difference_means_prev <- null_difference_means_current
    
  }
  return(list(pvalues, null_difference_means_current))
}


calcNumSimsNecessary <- function(evalAt, simulations){
  x <- simulations[1:evalAt]
  alpha = .1
  x0 = mean(x)
  s0 = sd(x)
  simsNecessary = (1.96 * s0 / x0 * alpha)^2 #246 of these 10
  return(simsNecessary)
}

#for(i in seq(1:ncol(null_difference_means))){
#  simsNecessary <- sapply(seq(1000, 10000, 1000), calcNumSimsNecessary, null_difference_means[,i])
#  plot(seq(1:length(simsNecessary)), simsNecessary, main = paste("Sims necessary for taxa ", i), xlab = "Number simulations done in 1000s", 
#       ylab = "Number of times you should do that number of simulations to get appropriate power")
#}

