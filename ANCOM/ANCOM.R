#Running ANCOM
setwd("~/Lab/Microbiome_ASD_16S/")
library(phyloseq)
install.packages("doParallel") 
install.packages("DT") 
install.packages("exactRankTests") 
install.packages("foreach")
install.packages("ggplot2") 
install.packages(“Rcpp”)
install.packages("openxlsx")
install.packages("shiny")
library(ape)
library("ancom.R")


#cleaning up/filtering
biom <- import_biom("data/otu_table_qc_age.biom")
otu <- otu_table(biom, taxa_are_rows = T)
map <- read.delim("mapping/mapping_final_withOutliers.txt")
#for the normalization later 
map$Treatment<-relevel(map$Treatment, ref = "Control")
rownames(map) <- map$SampleID
taxa<- read.csv("dada2/taxa.txt", row.names=1)
tree <- read.tree("tree/tree_mle_gtr_rootArchae_127.tre")
ps <- phyloseq(otu, sample_data(map), tax_table(as.matrix(taxa)), tree)
duplicates <- c("178.6", "179.6", "185.6", "188.6", "220.6", "221.6", "385") # choose to keep 385.6 b/c more reads than its duplicate
ps <- subset_samples(ps, !(SampleID %in% duplicates))
veryLowDiversity <- c("173", "367") #Outliers
ps <- subset_samples(ps, !(SampleID %in% veryLowDiversity))
phy_tree(ps) <- read.tree("tree/tree_mle_gtr_rootArchae_127.tre")


###Declare function to normalize
deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ Treatment)
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table(ps), phy_tree(ps))
  return(ps_deSeq)
}

###Declare function to filter
#```{r func_filterPrevalence}
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}
#```

###Declare function to make boxplots

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
    geom_jitter(aes(x = Group, y = value), position=position_jitter(0.2), cex=1.5, color="gray44") + 
    facet_wrap(~ variable  + pval, scale = "free") + labs(title = title, x = xlab, y = ylab) + scale_y_log10() + theme_minimal()
  
  print(p)
}

#Let's try with no filer 
#res_table<-otu_table(ps)
#res_table<-as.data.frame(res_table)
#res_table<-t(res_table)
#res_table<-as.data.frame(res_table) #eeds to be there twice
#metada_ps<-sample_data(ps)
#metada_ps<-as.data.frame(metada_ps)
#metada_ps$Treatment<-as.character(metada_ps$Treatment) 
#res_table$metada<-metada_ps$Treatment #same rown order so I can do that
#ancom.out.015 <- ANCOM( res_table, sig = 0.15, multcorr = 2, repeated=FALSE )

#checking results
#res_tax<-tax_table(ps)
#res_sign<-res_tax[rownames(res_tax) %in% ancom.out$detected]

#Now filter the prevalence
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
res_table_filt<-otu_table(ps_filt)
res_table_filt<-as.data.frame(res_table_filt)
res_table_filt<-t(res_table_filt)
res_table_filt<-as.data.frame(res_table_filt)
metada_ps<-sample_data(ps)
metada_ps<-as.data.frame(metada_ps)
res_table_filt$metada<-metada_ps$Treatment[rownames(metada_ps)%in%rownames(res_table_filt)]#same rown order so I can do that
ancom.out.filt.015 <- ANCOM(res_table_filt, sig = 0.15, multcorr = 2, repeated=FALSE )
ancom.out.filt.0.05 <- ANCOM(res_table_filt, sig = 0.05, multcorr = 2, repeated=FALSE )
saveRDS(ancom.out.filt.015, file="ANCOM/ancom.out.filt.015.rds")
saveRDS(ancom.out.filt.0.05, file="ANCOM/ancom.out.filt.0.05.rds")

#only two of them 
#plot the data after normalization 
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
ps_norm<-deSeqNorm(ps_filt)
esv_table_filt_norm<-otu_table(ps_norm)
esv_table_filt_norm<-esv_table_filt_norm[rownames(esv_table_filt_norm) %in% ancom.out.filt.015$detected,]

#boxplotting 
make_boxplots(t(esv_table_filt_norm), metada_ps$Treatment,c(0.15,0.15))

#which family are they from?
taxa_ps_norm<-tax_table(ps_norm)
tmp<-taxa_ps_norm[rownames(taxa_ps_norm) %in% ancom.out.filt.015$detected]
as.character(tmp[1,])
as.character(tmp[2,])
