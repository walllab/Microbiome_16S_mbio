#packags
library(phyloseq)
require(scales)
library(gridExtra)
library(ggplot2)
library(reshape2)

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

#ggtree(MyTree_achea)
ps_filt_norm<-readRDS("data/ps_filt_norm_117.RDS")
final_results<-read.csv("Summary_biomarkers/summ_results.csv", stringsAsFactors = F, row.names = 1)
MyTree_achea <- phy_tree(ps_filt_norm)
res.ps<-prune_taxa(rownames(final_results),ps_filt_norm)

#removing the taxa where we could not assign Genus: 
tax_res.ps<-as.data.frame(tax_table(res.ps))
tax.res.ps_genus<-tax_res.ps[!tax_res.ps$Genus %in% "g__" & !is.na(tax_res.ps$Genus),]
res.ps_genus<-prune_taxa(rownames(tax_res.ps),res.ps)
res.ps_genus_1<-transform_sample_counts(res.ps_genus,function(x) (x+1))

one_over_trans<- function() trans_new("one_over", function(x) log10(x+1), function(x) 10*exp(x)-1)

#First graph for a global look that we did not use in the paper
pdf("Summary_biomarkers/Biomarkers_abundance.pdf",width=10,height=6)
p1<-plot_bar(res.ps_genus_1, fill="Order", facet_grid=~Treatment, "Genus") 
p2<- p1 + coord_trans(y="one_over") + scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000))
p2 + scale_fill_manual(values=c("deeppink","darkorchid1","darkseagreen1","darkorange","gold1","lightskyblue", "blue","green"))
p2 + boxplot()
dev.off()    

#Ok let's try to do something easier to read for the paper: 
make_boxplots2 <- function(rsv_table, grouping, pvals, title = "", xlab = "", ylab = ""){
  pvals[pvals <= .01] = "***"
  pvals[pvals <= .05 & pvals > .01] = "**"
  df_grouped <- data.frame(cbind(rsv_table), grouping) # Create dataframe with RSV abundances and Group
  colnames(df_grouped) <- c(colnames(rsv_table), "Group") # Rename the columns so we can refer to them later
  grouped_stacked <- melt(df_grouped, id = "Group") # Put dataframe in form that's easy to use with ggplot
  rsv_table[rsv_table == 0] = .00001 #adding a small value for the plot with log transform 
  # Include Genus name in dataframe for graph labelling
  match_seq_to_pval <- data.frame(rsv_name = colnames(rsv_table), pval_adj = pvals) # Create little mapping dataframe for rsv_names to their genuses
  grouped_stacked$pval <- as.character(match_seq_to_pval$pval_adj[match(grouped_stacked$variable, match_seq_to_pval$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  # Plot! The function facet_wrap will break graphs up by whatever variables you put after the '~' character. In this case, we want to break it up by RSV name AND genus name
  p.list <- lapply(sort(unique(grouped_stacked$variable)), function(i) {
      ggplot(grouped_stacked[grouped_stacked$variable==i,], aes(x=Group, y = value)) + 
      geom_boxplot(aes(fill = Group)) +
      geom_jitter(aes(x = Group, y = value, color= factor(Group)), position=position_jitter(0.2), cex=1.5, alpha=.3, colour="grey44") + 
      #facet_wrap(~ variable  + pval, scale = "free") + 
      labs(title = i, y = "Log10 counts normalized by Size Factors", x="") + 
      scale_y_log10() + 
      #theme_minimal(base_size=2) +
      scale_fill_manual(values=alpha(c("red","blue"),0.3), drop = FALSE) +
      #guides(fill=FALSE) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=1),plot.title = element_text(size = 12) )
  })
  return(p.list)
}

#now running the function above
tmp_table <- t(otu_table(res.ps))
tmp_table<-as.data.frame(tmp_table)
#make sure that the final results is on the same order
final_results <- final_results[colnames(tmp_table), ]
#final_results_maude <- final_results[order(colnames(tmp_table)),] #ordered incorrectly
colnames(tmp_table) <- paste(final_results$ESV,final_results$Genus)

#remove NA 
colnames(tmp_table) <-gsub("NA","",colnames(tmp_table))
#and now order for the graph 
order_for_graph<-sapply(colnames(tmp_table), function(x) substr(x,5,6))
order_for_graph<-gsub(" ","",order_for_graph)
order_for_graph<-as.numeric(order_for_graph)
tmp_table<-tmp_table[,order(order_for_graph)]

#make SURE that the grouping of the function is in the right order
bla<-as.data.frame(sample_data(res.ps))
bla<-bla[order(rownames(tmp_table)),]

p.list<-make_boxplots2(tmp_table, bla$Treatment, final_results$pvals_adj)

#rsv_table<-tmp_table
#grouping<-sample_data(res.ps)$Treatment
#pvals<-final_results$pvals_adj

pdf("Summary_biomarkers/plot_each_taxa.pdf", width=20, height=15)

lay = rbind(c(1,2,3,11,12,13), 
            c(4,5,6,14,15,16),
            c(7,8,9,17,18,19),
            c(10,22,23,20,21,24))
grid.arrange(grobs = p.list , layout_matrix = lay)

dev.off()
#Still having a problem with the colors when the first group has no data: thent he red shift to right.
#I'm fixing this manualy with inkscape  

for (i in 1:length(rownames(final_results))){
  tax_table(ps_filt_norm)[rownames(tax_table(ps_filt_norm)) == rownames(final_results)[i], 4:6]
}
final_results$Genus[i]

