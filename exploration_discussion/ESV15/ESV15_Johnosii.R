library(phyloseq)  
library(ape)
library(dplyr)
library(reshape2)
library(ggplot2)
library(igraph)
require(scales)
library(gridExtra)
library(bitops)
library(RCurl)

setwd("~/Lab/Microbiome_ASD_16S/")

############################################################
############################################################
#ESV : uniformis analysis of specifics to this strain
#extracion of the genoms mathincg 100% rRna from ESV , and not the closest taxa, using IMG blast
############################################################
#check is there are unique ESV here
############################################################

same_genomes<- function(path_closets_metabo_csv,path_ESV_metolite_csv){
  metabolites_pathways_closest_ESV<-read.csv(path_closets_metabo_csv, row.names = 1, stringsAsFactors = F)
  metabolites_pawthays_ESV<-read.csv(path_ESV_metolite_csv, row.names = 1, stringsAsFactors = F)
  tmp<-length(setdiff(colnames(metabolites_pathways_closest_ESV),colnames(metabolites_pawthays_ESV)))
  tmp_2<-length(setdiff(colnames(metabolites_pawthays_ESV),colnames(metabolites_pathways_closest_ESV)))
  return(min(tmp,tmp_2))
}


creating_df_with_both<- function(path_closets_metabo_csv,path_ESV_metolite_csv){
  metabolites_pathways_closest_ESV<-read.csv(path_closets_metabo_csv, row.names = 1, stringsAsFactors = F)
  colnames(metabolites_pathways_closest_ESV)[2:dim(metabolites_pathways_closest_ESV)[2]]<-gsub("^","CLOS_",colnames(metabolites_pathways_closest_ESV)[2:dim(metabolites_pathways_closest_ESV)[2]])
  metabolites_pawthays_ESV<-read.csv(path_ESV_metolite_csv, row.names = 1, stringsAsFactors = F)
  #remove the names in one of th two since we use KO for cbind
  metabolites_pawthays_ESV<-metabolites_pawthays_ESV[colnames(metabolites_pawthays_ESV) != "Name"] 
  ESV_study<-cbind(metabolites_pathways_closest_ESV,metabolites_pawthays_ESV)
  ESV_study<-ESV_study[rowSums(ESV_study[,2:length(ESV_study)])!=0,] #453  KOs
  #saveing the names and removing it form the df 
  ESV_study_names<-as.data.frame(ESV_study$Name, stringsAsFactors=F)
  rownames(ESV_study_names)<-rownames(ESV_study)
  ESV_study<-ESV_study[,-1]
  to_return<-list()
  to_return[[1]]<-ESV_study_names
  to_return[[2]]<-ESV_study
  return(to_return)
}


#Is there any KO that is present in ESV and not in the rest 

finding_ko_present<-function(ESV_study,ESV_study_names){
  #replace the Ko by 1 or zero 
  ESV_study[ESV_study>0]<-1
  if (is.null(dim(ESV_study[,grep("CLOS_", colnames(ESV_study))]))){
    ESV_only_present<-rownames(ESV_study)[ESV_study[,grep("CLOS_", colnames(ESV_study))] == 0]
  } else {
    ESV_only_present<-rownames(ESV_study)[rowSums(ESV_study[,grep("CLOS_", colnames(ESV_study))]) == 0]
  }
  if (is.null(dim(ESV_study[,grep("CLOS_", colnames(ESV_study), invert =T)]))){
    ESV_only_not_present<-rownames(ESV_study)[ESV_study[,grep("CLOS_", colnames(ESV_study), invert =T)] == 0]
  } else {
    ESV_only_not_present<-rownames(ESV_study)[rowSums(ESV_study[,grep("CLOS_", colnames(ESV_study), invert =T)]) == 0]
  }
  #what are they? 
  Summary_df<-as.data.frame(c(ESV_only_present,ESV_only_not_present), stringsAsFactors=F)
  tmp_a<-rep("pres",length(ESV_only_present))
  tmp_b<-rep("abs",length(ESV_only_not_present))
  Summary_df$pres_abs<-c(tmp_a,tmp_b)
  colnames(Summary_df)<-c("KO","PRES_ABS")

  #adding the number of strain with it 
  if (is.null(dim(ESV_study[Summary_df$KO[Summary_df$PRES_ABS == "pres"],grep("CLOS_", colnames(ESV_study), invert=T)]))){
    pres_number<-ESV_study[Summary_df$KO[Summary_df$PRES_ABS == "pres"],grep("CLOS_", colnames(ESV_study), invert=T)]
  } else {
    pres_number<-rowSums(ESV_study[Summary_df$KO[Summary_df$PRES_ABS == "pres"],grep("CLOS_", colnames(ESV_study), invert=T)])
  }
  if (is.null(dim(ESV_study[Summary_df$KO[Summary_df$PRES_ABS == "abs"],grep("CLOS_", colnames(ESV_study))]))){
    abs_number<-ESV_study[Summary_df$KO[Summary_df$PRES_ABS == "abs"],grep("CLOS_", colnames(ESV_study))]
  } else {
    abs_number<-rowSums(ESV_study[Summary_df$KO[Summary_df$PRES_ABS == "abs"],grep("CLOS_", colnames(ESV_study))])
  }

  
  Summary_df$number_genome_present<-c(pres_number,abs_number)
  #and finally the total number of genome consider in each case
  tmp_a<-rep(length(grep("CLOS_", colnames(ESV_study), invert=T)),length(ESV_only_present))
  tmp_b<-rep(length(grep("CLOS_", colnames(ESV_study))),length(ESV_only_not_present))
  Summary_df$nb_genome_considered<-c(tmp_a,tmp_b)
  #fixing the rownames
  #Now adding the name 
  Summary_df<-cbind(ESV_study_names[Summary_df$KO,], Summary_df)
  #final little edit
  Summary_df$KO<-gsub("KO:","",Summary_df$KO)
  return(Summary_df)}


#and the API
run_API_KEGG<-function(Summary_df){
  #working on the API
  pathlist_sub=list()
  pathlist_prod=list()
  to_return=list()
  for (i in 1:length(Summary_df$KO)){
    cat(i, "\t",Summary_df$KO[i],"\t")
    reac<-c()
    reac <-getURL(paste("http://rest.kegg.jp/link/reaction/",Summary_df$KO[i],sep=""))
    cat(i, "\t",reac,"\t")
    reac<-strsplit(reac,"\n")
    reac <-strsplit(reac[[1]],"\t")
    
    sub_tmp<-c()
    produc_tmp<-c()
    for (j in (1:length(reac))){
      reac1<-c()
      reac1 <-grep("rn:",reac[[j]], value = T)
      reac1<- gsub("rn:","",reac1)
      reac1<- gsub("\n","",reac1)
      comp_tmp1<-c()
      comp_tmp1<-getURL (paste("http://rest.kegg.jp/list/",reac1, sep=""))
      comp_tmp1<-strsplit(comp_tmp1,";")
      comp_tmp1<-strsplit(comp_tmp1[[1]][2], "<=>")
      sub_tmp1<-comp_tmp1[[1]][1]
      sub_tmp<-c(sub_tmp,sub_tmp1)
      cat(sub_tmp,"\n")
      produc_tmp1<-gsub("\n","",comp_tmp1[[1]][2])
      produc_tmp<-c(produc_tmp,produc_tmp1)
    }
    pathlist_prod[[i]]<-produc_tmp
    pathlist_sub[[i]]<-sub_tmp
  }
  pathlist_prod<-sapply(pathlist_prod, function(x){paste(x,collapse = " ;")})
  pathlist_sub<-sapply(pathlist_sub, function(x){paste(x,collapse = " ;")})
  to_return[[1]]<-pathlist_prod
  to_return[[2]]<-pathlist_sub
  return(to_return)}

#####################################################################################################################################################################
#####################################################################################################################################################################
ps <- readRDS("data/ps_filt_norm_117_BS100.RDS")
ps.clos<-subset_taxa(ps, Family=="f__Moraxellaceae") #cn't run it: only one from this family! 


#then find the two that we have for ESV
plot_tree(ps.clos)
df_rel<-as.data.frame(tax_table(ps.clos))
####Have a loog at our biomarkers 
final_results<-read.csv("Summary_biomarkers/summ_results.csv", stringsAsFactors = F, row.names=1)
rownames(for_plot)<-for_plot$Column1
for_plot<-merge(final_results, df_rel,by =0, all.y=T)
for_plot$enriched<-as.factor(for_plot$enriched)
#removing the column of the taxa not needed 
for_plot<-for_plot[,-c(2:8)]
rownames(for_plot)<-for_plot$Row.names
for_plot<-as.matrix(for_plot)
tax_table(ps.clos)<-for_plot

#first find the genomes in IMG that exactly matches it: 
tree_analysis<-plot_tree(ps.clos,color="enriched", label.tips = "Genus.y", title = "Family Lachnospiraceae")

node_of_interets<-tree_analysis$data$V1[tree_analysis$data$OTU == rownames(final_results)[final_results$ESV == "ESV#15"]]
node_of_interets<-node_of_interets[!is.na(node_of_interets)]
tmp<-tree_analysis$data$OTU [tree_analysis$data$V2[tree_analysis$data$V1 == node_of_interets]]
to_extract_with_ESV <-tmp[tmp != rownames(final_results)[final_results$ESV == "ESV#15"]]


###Went into IMG and generated the file "closest_ESV_KEGG_metabolites.csv" and "exploration_discussion/ESV_KEGG_metabolites.csv"
#Running this for ESV
#STEP1
same_genomes("exploration_discussion/ESV15/closest_ESV15_metabolites.csv","exploration_discussion/ESV15/ESV15_metabolites.csv") #check that is different from wero
ESV_to_parse<-creating_df_with_both("exploration_discussion/ESV15/closest_ESV15_metabolites.csv","exploration_discussion/ESV15/ESV15_metabolites.csv")
#Now function to analyse and extract the ones that are unique 
ESV_study<-ESV_to_parse[[2]]
ESV_study_names<-ESV_to_parse[[1]]

#STEP 2
Summary_df<-finding_ko_present(ESV_study,ESV_study_names)

#STEP 3: API
#And now adding it to the dataframe
KEGG_API_ESV<-run_API_KEGG(Summary_df)
#Summary_df$sub<-pathlist_sub
#Summary_df$prod<-pathlist_prod
Summary_df$sub<-KEGG_API_ESV[[2]]
Summary_df$prod<-KEGG_API_ESV[[1]]

#And saving it all
write.csv(Summary_df, file="exploration_discussion/ESV15/ESV15_specificity.csv")
write.csv(ESV_study_names, file="exploration_discussion/ESV15/ESV15_study_names.csv")
write.csv(ESV_study, file="exploration_discussion/ESV15/ESV15_study.csv")

##############################################################################################################
##############################################################################################################
#Tree for the paper, more information than just one branch for the visual 
ps.clos<-subset_taxa(ps, Family =="f__Moraxellaceae")

#then find the two that we have for ESV
plot_tree(ps.clos)
df_rel<-as.data.frame(tax_table(ps.clos))
####Have a loog at our biomarkers 
final_results<-read.csv("Summary_biomarkers/summ_results.csv", stringsAsFactors = F, row.names=1)
rownames(for_plot)<-for_plot$Column1
for_plot<-merge(final_results, df_rel,by =0, all.y=T)
for_plot$enriched<-as.factor(for_plot$enriched)
#removing the column of the taxa not needed 
for_plot<-for_plot[,-c(2:8)]
rownames(for_plot)<-for_plot$Row.names
for_plot<-as.data.frame(for_plot)
tmp<-apply(for_plot,1, function(x) {paste(x[names(x) == "Genus.y"],x[names(x) == "ESV"],collapse="_")})
for_plot$names<-tmp
for_plot$names<-gsub("NA","",tmp)
for_plot<-as.matrix(for_plot)
tax_table(ps.clos)<-for_plot
plot_tree(ps.clos,color="enriched", label.tips ="names" , title = "Family Lachnospiraceae")


