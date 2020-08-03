library(phyloseq)  
library(ape)
library(dplyr)
library(structSSI)
library(reshape2)
library(ggplot2)
library(igraph)
library(readr)
library(gage)

setwd("~/Lab/Microbiome_ASD_16S/")
set.seed(100)
source('~/Lab/Microbiome_ASD_16S/pair_analysis/mc.R')
source('~/Lab/Microbiome_ASD_16S/scripts/getTables.R')
source('~/Lab/Microbiome_ASD_16S/metagenomseq/ZIG_mixture_model.R')
source('~/Lab/Microbiome_ASD_16S/pathway_enrichment/pval_plotting.R')

###Load Data
ps <- readRDS("data/ps_noDuplicates.RDS")

phy_tree(ps) <- read.tree("tree/tree_mle_gtr_rootArchae_127.tre")

#################################################Declaring all of our functions#######################
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

#Remove the samples which don't have ASD and their sibling according to video classifier and MARA:
ID_not_ASD<-sample_data(ps)$SampleID[sample_data(ps)$classifier > 0 & sample_data(ps)$Treatment == "Aut"]
#and find thier siblings 
pair_sibling_rem<-sample_data(ps)$Pair[sample_data(ps)$SampleID %in% ID_not_ASD]
samples_to_remove_not_ASD<-sample_data(ps)$SampleID[sample_data(ps)$Pair %in% pair_sibling_rem]
sample_to_keep<-rownames(sample_data(ps))[!rownames(sample_data(ps)) %in% samples_to_remove_not_ASD]
ps <-prune_samples(sample_to_keep,ps)
###Filter and normalize
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
ps_filt_norm <- deSeqNorm(ps_filt)

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

################################################# Pair Analysis #######################

###Run permutation test
ps_pairs <- subset_samples(ps_filt_norm, !(SampleID %in% c("45.1", "233"))) #take out the third siblings
pairs_to_keep <- sample_data(ps_pairs)$Pair[duplicated(sample_data(ps_pairs)$Pair)]
ps_pairs <- subset_samples(ps_pairs, Pair %in% pairs_to_keep)

#Pick whether to load null distribution, or re-run
#nullDists <- readRDS(paste("pair_analysis/nullDist_", prevFiltThresh, "filter_", nsamples(ps_pairs), "s_adj.RDS", sep = ""))
#nullDists <- runSimulations(otu_table(ps_pairs), sample_data(ps_pairs), numSims= 10000) #returns a simulation X factor table
#saveRDS(nullDists, paste("pair_analysis/nullDist_", prevFiltThresh, "filter_", nsamples(ps_pairs), "s_adj.RDS", sep = ""))
nullDists<-readRDS("pair_analysis/nullDist_0.03filter_110s_adj.RDS")
mc_res = mcPermutationTest(ps_pairs, null_difference_means = nullDists, numSims = 10000) #Pass loaded null dists
mc_res$nb_non_zero = apply(otu_table(ps_pairs)[mc_res$seqs, ], 1, function(tax) return(sum(as.numeric(tax)>0))) #find the number of samples that taxa is present in
saveRDS(mc_res, file= paste("pair_analysis/mc_", prevFiltThresh, "filter_", nsamples(ps_pairs), "s_adj.RDS", sep = ""))

#Visualize with boxplots
tmp_table <- t(otu_table(ps_pairs)[mc_res$seqs, ])
colnames(tmp_table) <- paste(mc_res$Family, mc_res$Genus, mc_res$Species)

make_boxplots(tmp_table, sample_data(ps_pairs)$Treatment, mc_res$pvals_adj)
rsv_table<- tmp_table
grouping<-sample_data(ps_pairs)$Treatment
pvals<-mc_res$pvals_adj


#Check out distributions: how believable is that p value?
pdf("pair_analysis/mc_distributions.pdf",width=6,height=4,paper='special') 
actual_means = getDifferenceMeans(otu_table(ps_pairs), sample_data(ps_pairs))
rownames(mc_res) <- mc_res$seqs
for(taxa in mc_res$seqs){
  #Get actual line and pvalue
  pval <- mc_res[taxa, ]$pvals_adj
  actual_mean <- actual_means[taxa]
  
  #Get taxa name
  taxName <- paste(mc_res[taxa, ]$Family, mc_res[taxa, ]$Genus, mc_res[taxa, ]$Species, collapse = " ")
  
  #Graph 
  numSims = 10000
  hist(nullDists[1:numSims, taxa], prob = F, main = paste("Null Distribution obtained by ", numSims, "Permutations \n", "Taxa: ", taxName), 
       xlab = "Value of sibling difference mean",
       ylab = "Frequency", breaks = 10)
  #lines(density(nullDists[1:numSims, taxa]), col = "red")
  abline(v = actual_mean, col = "yellow")
  
}
dev.off()

#################################################Run Deseq#######################

###Run DESeq proper
runDESeq <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ Treatment) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("Treatment", "Aut", "Control"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
  }

deseq_res <- runDESeq(ps_filt)
deseq_res$enriched <- ifelse(deseq_res$log2FoldChange > 0,"Aut", "Control")
saveRDS(deseq_res, file="DESeq/deseq_res_117.rds")


#make_boxplots(t(otu_table(ps_filt_norm)[rownames(deseq_res), ]), sample_data(ps_filt_norm)$Treatment, paste(deseq_res$Family, deseq_res$Genus, deseq_res$Species), deseq_res$padj)

#################################################ZIG with metagenomSeq#######################

###Run ZIG model fitting and prediction
zig_res <- run_metagenom_seq(ps_filt,30)
zig_res <- data.frame(cbind(zig_res), tax_table(ps_filt)[rownames(zig_res), ])
zig_res$enriched <- ifelse(zig_res$TreatmentControl < 0, "Aut", "Control")
saveRDS(zig_res, file="metagenomseq/zig_res_117.rds")

#################################################Piphillin analysis#######################

#Run this to redo calculations. Here, we load the saved object to save time
#loc.koTxt <- "pathway_enrichment/piphillan_results_noDuplicates_filter03_deseq/kegg_output/ko_abund_table_unnorm.txt"
#table <- as.data.frame(read_delim(loc.koTxt, "\t", escape_double = FALSE, trim_ws = TRUE))
#rownames(table) <- table[,1]
#table <- table[,-1]

#create pathway table by summing the relative contribution of each KO to each pathway
#Divide each KO abundance by the number of pathways that it is a part of, based on the database
#Don't need to re-run, just load the saved one
#numPathways_perKO <- sapply(rownames(table), function(ko_name) return(sum(unlist(gs) %in% ko_name)))
#numPathways_perKO <- numPathways_perKO[numPathways_perKO != 0]
#table <- table[rownames(table) %in% names(numPathways_perKO),]
#table_norm <- table / numPathways_perKO
#saveRDS(table_norm, "pathway_enrichment/ko_abund_table_normByPathway.Rds")

#Run and plot GSEA for pathways
gs <- readRDS("pathway_enrichment/pathways_gsea.rds")
ps_piph <- readRDS("data/ps_noDuplicates.RDS")

table_norm <- readRDS("pathway_enrichment/ko_abund_table_normByPathway.Rds")
colnames(table_norm) <- sample_data(ps_piph)$SampleID

#Remove samples that don't have autism and their siblings
sample_to_keep <- rownames(sample_data(ps_piph))[!rownames(sample_data(ps_piph)) %in% samples_to_remove_not_ASD]
table_norm <- table_norm[ , sample_to_keep]



ps_kos <- phyloseq(otu_table(table_norm, taxa_are_rows = T), sample_data(ps_filt_norm))
ps_kos_aut <- subset_samples(ps_kos, Treatment == "Aut")
ps_kos_control <- subset_samples(ps_kos, Treatment == "Control")

exp_table <- cbind(otu_table(ps_kos_control), otu_table(ps_kos_aut))
map <- rbind(sample_data(ps_kos_control), sample_data(ps_kos_aut))

e <- gage(exp_table, gsets = gs, 
          ref = seq(1,nsamples(ps_kos_control)), 
          samp = seq((nsamples(ps_kos_control) + 1),(nsamples(ps_kos_control)+ nsamples(ps_kos_aut))),
          compare = 'unpaired')

greater<-as.data.frame(e$greater)[as.data.frame(e$greater)$q.val<0.15 & !is.na(as.data.frame(e$greater)$q.val),]
less<-as.data.frame(e$less)[as.data.frame(e$less)$q.val<0.15 & !is.na(as.data.frame(e$less)$q.val),]

#save the results 

relevant_GSEA<-rbind(less[order(less$q.val),1:5],greater[order(greater$q.val),1:5])
saveRDS(relevant_GSEA,file="pathway_enrichment/relevant_GSEA_110.rds")
pdf("pathway_enrichment/GSEA_relevant.pdf")
plotPValues(e$greater, e$less, condition="All Autism and Control \n KO pathways")
dev.off()


#Put results together 
zig_res_117<-readRDS("metagenomseq/zig_res_117.rds")
mc_03filter_110s_deseq_adj<-readRDS("pair_analysis/mc_0.03filter_110s_adj.RDS")
deseq_res<-readRDS("DESeq/deseq_res_117.rds")
write.csv(zig_res,file="comparing_results/zig_res.csv")
write.csv(mc_03filter_110s_deseq_adj,file="comparing_results/mc_03filter_110s_deseq_adj.csv")
write.csv(deseq_res,file="comparing_results/deseq_res.csv")


####ML####
library(e1071)
library(randomForest)
library(caret)
library(pROC)
set.seed(10)
ml <- function(ps, model = "randomForest", seqs){
  abund <- otu_table(ps)
  accs <- c()
  flds <- createFolds(unique(sample_data(ps)$Pair), k = 10, list = TRUE, returnTrain = FALSE)
  for(i in seq(1, length(flds))){
    training_pairs <- as.numeric(unlist(flds[-i]))
    testing_pairs <- flds[[i]]
    training_bools <- sample_data(ps)$Pair %in% training_pairs
    testing_bools <- sample_data(ps)$Pair %in% testing_pairs
    
    if(model == "randomForest"){
      trained_model <- randomForest(t(abund[seqs, training_bools]), y = as.factor(sample_data(ps)$Treatment[training_bools]), type = "classification")
    }
    if(model == "svm"){
      trained_model <- svm(t(abund[seqs, training_bools]), y = as.factor(sample_data(ps)$Treatment[training_bools]), type = "C-classification")
    }
    
    preds <- predict(trained_model, t(abund[seqs, testing_bools]))
    print(preds)
    acc <- sum(preds == as.factor(sample_data(ps)$Treatment[testing_bools])) / length(preds)
    accs <- c(accs, acc)
    ROC <- roc(factor(sample_data(ps)$Treatment[testing_bools], levels = c("Control", "Aut")), factor(preds, levels = c("Control", "Aut"), ordered = T))
    if(model == "randomForest"){
      title = "Random Forest 10 fold Cross Validation ROC curves, unsmoothed"
    }
    if(model == "svm"){
      title = "SVM 10 fold Cross Validation ROC curves, unsmoothed"
    }
    if(i == 1){
      plot.roc(ROC, add = F, main = title)
    }else{
      plot.roc(ROC, add = T, main = title)
    }
    
  }
  print(mean(accs)) #62% : better than random!
  
}

ml(ps_filt_norm, model = "svm", seqs = c(as.character(mc_res$seqs), rownames(zig_res), rownames(deseq_res))) #64% accuracy
ml(ps_filt_norm, model = "randomForest", seqs = c(as.character(mc_res$seqs), rownames(zig_res), rownames(deseq_res)))

#We don't get good accuracy due to the heterogeneity of autism: Realistically, we should be able to subset and get much better results, however, our sample size 
#prohibits us from doing this! Hence, phase 2!
