#We want to see where each of the results actually plot on the tree 
library(ggtree)
library(ggplot2)
library(devtools)
require(gridExtra)
library(grid)
library(lattice)
library(cowplot)
library(gtable)
library(ddplyr)
library(scales)
#source("https://bioconductor.org/biocLite.R") #Some problem between phyloseq and lastest ape version solved in phyloseq_1.22.3.tgz
#biocLite("phyloseq")
library(phyloseq)

setwd("~/Lab/Microbiome_ASD_16S/Microbiome_ASD_16S")
#not needed as already in the phyloseq object
#MyTree_achea <- read.tree("tree/tree_mle_gtr_rootArchae_127_BS100.tre")
#phy_tree(ps) <-MyTree_achea
final_results<-read.csv("Summary_biomarkers/summ_results.csv", stringsAsFactors = F)
#results are consistent between the two tress, we keep the archea tree
ps<-readRDS("data/ps_filt_norm_117_BS100.RDS")

#the taxa_to_plot

marg<-c(0.1,0.19,0.1,0.1,0.2,0.3,0.1,0.1,0.1,0.1)
spac<-c(0.01,0.009,0.1,0.1,0.05,0.1,0.01,0.01,0.1,0.1)

order_to_map<-c("o__Bacteroidales","o__Clostridiales","o__Pasteurellales","o__Coriobacteriales", "o__Desulfovibrionales", "o__Bifidobacteriales", "o__Erysipelotrichales")
title_trees<-c("Bacteroidales","Clostridiales","Pasteurellales","Coriobacteriales", "Desulfovibrionales", "Bifidobacteriales", "Erysipelotrichales")
#"o__Pseudomonadales" only one ESV in this order: no need to plot the tree 
#

#we need specific marging here
plot_order_mapped=list()
for (i in 1:length(order_to_map)){
  ps.clos <- subset_taxa(ps, Order==order_to_map[i])
  df_rel<-as.data.frame(tax_table(ps.clos))
  
  enriched_aut_bool <- rownames(df_rel) %in% final_results$Column1[final_results$enriched == "Aut"]
  enriched_nt_bool <- rownames(df_rel) %in% final_results$Column1[final_results$enriched == "Control"]
  
  df_rel$Taxa<-rep("Non relevant",dim(df_rel)[1])
  df_rel$Taxa[enriched_aut_bool]<-"ASD biomarkers"
  df_rel$Taxa[enriched_nt_bool]<-"Neurotypical biomarkers"
  df_rel$Taxa<-as.factor(df_rel$Taxa)
  
  df_rel$to_annotate<-rep("",dim(df_rel)[1])
  
  ESV_aut <- final_results$ESV[match(rownames(df_rel)[enriched_aut_bool], final_results$Column1)]
  ESV_nt <- final_results$ESV[match(rownames(df_rel)[enriched_nt_bool], final_results$Column1)]
  
  df_rel$to_annotate[enriched_aut_bool]<- gsub("NA", "", paste(ESV_aut, 
                                                gsub("g__", "", as.character(df_rel$Genus[enriched_aut_bool])),
                                                gsub("s__", "", as.character(df_rel$Species[enriched_aut_bool])),
                                                sep=" "))
  df_rel$to_annotate[enriched_nt_bool]<- gsub("NA", "", paste(ESV_nt, 
                                              gsub("g__", "", as.character(df_rel$Genus[enriched_nt_bool])),
                                              gsub("s__", "", as.character(df_rel$Species[enriched_nt_bool])),
                                              sep=" "))
  df_rel$to_annotate<-as.factor(df_rel$to_annotate)
  
  df_rel$label_colors <- rep("Non relevant", nrow(df_rel))
  df_rel$label_colors[enriched_aut_bool] <- "Aut"
  df_rel$label_colors[enriched_nt_bool] <- "Control"
  tax_table(ps.clos)<- tax_table(as.matrix(df_rel))

plot_order_mapped[[i]] <- plot_tree_mod(ps.clos, ladderize = "left",
                                  color = "plot_colors", 
                                  title = order_to_map[i], 
                                  base.spacing = spac[i], 
                                  plot.margin = marg[i], 
                                  label.tips = "to_annotate",
                                  label.colors = "label_colors",
                                  text.size = 3, 
                                  nodelabf = nodeplotboot(highthresh=80, lowcthresh=50, size=2, hjust=1.8),
                                  size = "abundance") + 
  theme(legend.position="none", axis.line=element_blank(), plot.title = element_text(color="#666666", face="bold", size=12, hjust=0)) +                        
  scale_color_manual(values = c("Non relevant" = "gray", "Aut" = "firebrick1","Control"="dodgerblue1")) +
  scale_size_continuous(range = c(0.1, 1)) + 

  ggtitle(title_trees[i], subtitle = NULL)
} 


#and edo the last one with the legend
#for_legend<-plot_tree(ps.clos, 
#                      ladderize="left", 
#                      color="Taxa" , 
#                      title=order_to_map[3], 
#                      base.spacing =0.03,
#                      plot.margin= 1.5, 
#                      label.tips = "to_annotate", 
#                      text.size = 4) +
#  scale_color_manual(values = c("Non relevant" = "gray", "ASD biomarkers" = "firebrick1","Neurotypical biomarkers"="dodgerblue1"))

#plot_order_mapped[[8]]<-get_legend(for_legend)

lay <- rbind(c(2,2,1,1),
              c(2,2,1,1),
              c(2,2,3,3),
              c(2,2,4,4),
              c(2,2,7,7),
              c(2,2,6,6),
              c(2,2,5,5))

pdf("Summary_biomarkers/all_trees_treatment.pdf",width=12,height=10)
grid.arrange(grobs = plot_order_mapped, layout_matrix = lay)
dev.off()



