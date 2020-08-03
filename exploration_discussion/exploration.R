library(phyloseq)  
library(ape)
library(dplyr)
library(reshape2)
library(ggplot2)
library(igraph)
require(scales)
library(gridExtra)

setwd("~/Lab/Microbiome_ASD_16S/")

###Declare function to normalize
deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ Treatment )
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table(ps), phy_tree(ps))
  return(ps_deSeq)
}

###Declare function to filter
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

make_boxplots <- function(rsv_table, grouping, pvals, title = "", xlab = "", ylab = ""){
  # rsv_table should be a sample x feature table
  # grouping should match the samples and describe what condition the samples have
  # extraFeatureLabeling is a vector that matched the features in the rsv_table that describes extra info you want to display, in this case genus
  #pvals <- round(pvals, digits = 2)
  pvals[pvals < .01] = "***"
  df_grouped <- data.frame(cbind(rsv_table), grouping) # Create dataframe with RSV abundances and Group
  
  colnames(df_grouped) <- c(colnames(rsv_table), "Group") # Rename the columns so we can refer to them later
  
  grouped_stacked <- melt(df_grouped, id = "Group") # Put dataframe in form that's easy to use with ggplot
  
  # Include Genus name in dataframe for graph labelling
  #match_seq_to_extraInfo <- data.frame(rsv_name =colnames(rsv_table), extraInfo = extraFeatureLabeling) # Create little mapping dataframe for rsv_names to their genuses
  match_seq_to_pval <- data.frame(rsv_name = colnames(rsv_table), pval_adj = pvals) # Create little mapping dataframe for rsv_names to their genuses
  #grouped_stacked$extraInfo <- as.character(match_seq_to_extraInfo$extraInfo[match(grouped_stacked$variable, match_seq_to_extraInfo$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  grouped_stacked$pval <- as.character(match_seq_to_pval$pval_adj[match(grouped_stacked$variable, match_seq_to_pval$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  
  # Plot! The function facet_wrap will break graphs up by whatever variables you put after the '~' character. In this case, we want to break it up by RSV name AND genus name
  p <- ggplot(grouped_stacked, aes(x=Group, y = value)) + geom_boxplot(aes(fill = Group)) +
    geom_jitter(aes(x = Group, y = value), position=position_jitter(0.2), cex=1.5, color="gray44") + facet_wrap(~ variable  + pval, scale = "free") + labs(title = title, x = xlab, y = ylab) + scale_y_log10() + theme_minimal()
  
  print(p)
}

###Load Data
ps <- readRDS("data/ps_noDuplicates.RDS")
phy_tree(ps) <- read.tree("tree/tree_mle_gtr_rootArchae_127.tre")
###Filter and normalize
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
ps_filt_norm <- deSeqNorm(ps_filt)


############################################################
#Let's work with the ratio firmicute over Bacteroides
############################################################
taxa_df<-as.data.frame(tax_table(ps_filt_norm))
taxa_firmicutes<-taxa_df[taxa_df$Phylum=="p__Firmicutes",]
taxa_bacteroides<-taxa_df[taxa_df$Phylum=="p__Bacteroidetes",]
res.ps.firmicutes<-prune_taxa(rownames(taxa_firmicutes),ps_filt_norm)
res.ps.bacteroides<-prune_taxa(rownames(taxa_bacteroides),ps_filt_norm)
#and now the ratio 
ratio_formicute_over_bacteroides<-colSums(otu_table(res.ps.firmicutes))/colSums(otu_table(res.ps.bacteroides))
metada<-sample_data(res.ps.firmicutes)
wilcox.test(ratio_formicute_over_bacteroides[metada$Treatment=="Aut"], ratio_formicute_over_bacteroides[metada$Treatment=="Control"]) #no differences p-value = 0.5974, also barplot look really the same 
boxplot(ratio_formicute_over_bacteroides[metada$dairy == "None"], ratio_formicute_over_bacteroides[metada$dairy == "Dairy_int"]) #no differences p-value = 0.p-value = 0.6193, also barplot look really the same 


############################################################
#How about the diversity?
############################################################
#See diversity folder Chrisitne's analysis

