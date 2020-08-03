#Find the KOs from which genomes in the pathways 
setwd("Lab/Microbiome_ASD_16S/Microbiome_ASD_16S/")
gs <- readRDS("pathway_enrichment/pathways_gsea.rds") #all the patwhays from KEGG
ps_piph <- readRDS("data/ps_noDuplicates.RDS") #the taxa detected 
library(seqinr)

submitted_pi_phi<-read.fasta(file =("pathway_enrichment/repseqs.fasta"), as.string = TRUE)
submitted_pi_phi<-lapply(submitted_pi_phi, function(x) {x[1]})
submitted_pi_phi<-sapply(submitted_pi_phi, function(x) {toupper(x)})

#Now find the ESV of interest! 
final_results<-read.csv("Summary_biomarkers/summ_results.csv", stringsAsFactors = F, row.names=1)

submitted_pi_phi<-submitted_pi_phi[submitted_pi_phi %in% rownames(final_results)]

names(submitted_pi_phi) #gives you the otu with the KOs we need to grab 


ko_abund<-readRDS("pathway_enrichment/ko_abund_table_normByPathway.Rds")
load("pathway_enrichment/piphillan_output/450p_genome_contribution_table.RData")
otu_genome_hit_table<-read.table("pathway_enrichment/piphillan_output/200p_otu_genome_hit_table.txt", stringsAsFactors = F)
#find each genome with a hit for each ESV 
genomes_imp<-otu_genome_hit_table[otu_genome_hit_table$V1 %in% names(submitted_pi_phi),]
#we have only 6 ESVs that had genomes sequenced used in piphillin

#Which ESVs from biomarker list are able to be assigned to whole genome?
View(final_results[submitted_pi_phi[genomes_imp$V1],])



#Adding the sequence in another column 
tmp_b<-c()
for (i in 1:length(genomes_imp$V1)){
tmp<-submitted_pi_phi[names(submitted_pi_phi) == genomes_imp$V1[i]]
tmp_b<-c(tmp_b,tmp)
}
genomes_imp$V4<-tmp_b
#adding column name and saving the file 
colnames(genomes_imp)<-c("otu_piphilin","closest_KEGG_genome","%","ESV_Sequence")
#and finally the ESV names
ESV_names<-c("ESV1","ESV6","ESV21","ESV15","ESV14","ESV16")
genomes_imp$ESV_names<-ESV_names
#and maybe adding the 
write.csv(genomes_imp, file="pathway_enrichment/genomes_contribute_piph.csv")

#and now the KOs present in those genomes, save as a list: 
list_KO_per_genome<-list()
for (i in 1:length(genomes_imp$closest_KEGG_genome)){
list_KO_per_genome[[i]]<-genome_contribution$feature[genome_contribution$genome == genomes_imp$V2[i]]
}
names(list_KO_per_genome)<-genomes_imp$closest_KEGG_genome
#and now for each of the pathway
relevant_pathways<-readRDS("pathway_enrichment/relevant_GSEA_110.rds")
relevant_pathways<-relevant_pathways[relevant_pathways$q.val < 0.05,]

#And now the KOs involved in them:  
list_pathway_relevant_with_KOS<-gs[names(gs) %in% rownames(relevant_pathways)] #17 patwhays
#We need to have a list of KO with the genomes they belong to for each pathway with the genome
pathwas_with_KOs=list()
for (i in 1:length(list_pathway_relevant_with_KOS)){
  save_genome<-matrix(ncol=2,nrow=0)
  colnames(save_genome)<-c("KOs","genome")
  for (j in 1:length(list_KO_per_genome)){
    tmp_g<-list_pathway_relevant_with_KOS[[i]][list_pathway_relevant_with_KOS[[i]] %in% list_KO_per_genome[[j]]]
    tmp_g<-as.data.frame(tmp_g, stringsAsFactors=F)
    tmp_g$genome_name<-rep(names(list_KO_per_genome)[j],dim(tmp_g)[1])
    colnames(tmp_g)<-c("KOs","genome")
    save_genome<-rbind(save_genome,tmp_g)
  }
  pathwas_with_KOs[[i]]<-save_genome
}
names(pathwas_with_KOs)<-names(list_pathway_relevant_with_KOS)
###saveRDS(pathwas_with_KOs,"pathway_enrichment/pathwas_with_KOs.RDS")

#pathwas_with_KOs<-readRDS("pathway_enrichment/pathwas_with_KOs.RDS")
#Ok now let's explore the pathways of interest

#####################################################################
#Benzoate degradation ko00362
Benzoate_degradation_ko00362<-pathwas_with_KOs[grep("*ko00362*",names(pathwas_with_KOs))][[1]]
#some KOs are presetn in several 
Benzoate_degradation_ko00362_export<-as.data.frame(unique(Benzoate_degradation_ko00362$KOs), stringsAsFactors=F)
names_genome<-c()
for(i in 1:dim(Benzoate_degradation_ko00362_export)[1]){
  tmp<-Benzoate_degradation_ko00362$genome[Benzoate_degradation_ko00362$KOs == Benzoate_degradation_ko00362_export[i,1]]
  #adding the name of that genome
  tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
  names_genome<-c(names_genome,paste(tmp,collapse = "_"))
}
Benzoate_degradation_ko00362_export$genomes<-names_genome
#And now adding the color
Benzoate_degradation_ko00362_export$genomes<-as.factor(Benzoate_degradation_ko00362_export$genomes)
levels(Benzoate_degradation_ko00362_export$genomes) #looking ta it and determine the colors based on if they are biomarkers for which cohort
#[1] "ESV1_ESV15_ESV16"= "7733ff"(purple)  "ESV15"  = "33afff"(light blue)          "ESV16"  = 1130f4(medium blue)           "ESV21_ESV15_ESV14"="112392"(dark blue)  "ESV6" =f50808 (red)            
#[6] "ESV6_ESV15" == "7733ff"(purple)
levels(Benzoate_degradation_ko00362_export$genomes)<-c("#7733ff","#33afff","#1130f4","#112392","#f50808","#7733ff")
#Now exporting for kegg 
colnames(Benzoate_degradation_ko00362_export)<-c()
write.table(Benzoate_degradation_ko00362_export, file = "pathway_enrichment/Benzoate_degradation_ko00362/Benzoate_degradation_ko00362_export", quote = FALSE,sep = "\t", row.names =F )
#then exploring on KEGG GUI interface

##########################################################################################################################################
##########################################################################################################################################
#Butanoate metabolism ko00650
Butanoate_metabolism_ko00650<-pathwas_with_KOs[grep("*ko00650*",names(pathwas_with_KOs))][[1]]
Butanoate_metabolism_ko00650_export<-as.data.frame(unique(Butanoate_metabolism_ko00650$KOs), stringsAsFactors=F)
names_genome<-c()
KEGG_genomes_names<-c()
for(i in 1:dim(Butanoate_metabolism_ko00650_export)[1]){
  tmp<-Butanoate_metabolism_ko00650$genome[Butanoate_metabolism_ko00650$KOs == Butanoate_metabolism_ko00650_export[i,1]]
  #adding the name of that genome
  tmp_KEGG<-paste(tmp,collapse = "_")
  tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
  names_genome<-c(names_genome,paste(tmp,collapse = "_"))
  KEGG_genomes_names<-c(KEGG_genomes_names,tmp_KEGG)
}
Butanoate_metabolism_ko00650_export$genomes<-names_genome
Butanoate_metabolism_ko00650_export$genomes<-as.factor(Butanoate_metabolism_ko00650_export$genomes)
levels(Butanoate_metabolism_ko00650_export$genomes) 
Butanoate_metabolism_ko00650_export$genomes_save<-Butanoate_metabolism_ko00650_export$genomes

#we'll produce two maps, with ASD and without ASD
tmp_list<-gsub("ESV","",levels(Butanoate_metabolism_ko00650_export$genomes))
tmp_list<-strsplit(tmp_list,"_")
levels(Butanoate_metabolism_ko00650_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]<-rep("#b1b6bc",length(levels(Butanoate_metabolism_ko00650_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]))
#Now the ESV with 1-> 10
tmp_list<-gsub("ESV","",levels(Butanoate_metabolism_ko00650_export$genomes))
tmp_list<-as.numeric(tmp_list) #it's ok if there is a warning 
levels(Butanoate_metabolism_ko00650_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]<-rep("#f93a3a",length(levels(Butanoate_metabolism_ko00650_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]))
levels(Butanoate_metabolism_ko00650_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]<-rep("#529af0",length(levels(Butanoate_metabolism_ko00650_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]))

colnames(Butanoate_metabolism_ko00650_export)<-c()
write.table(Butanoate_metabolism_ko00650_export, file = "pathway_enrichment/Butanoate_metabolism_ko00650/Butanoate_metabolism_ko00650_export.tsv", quote = FALSE,sep = "\t", row.names =F )
#####################################################################################
#and now we will color per genomes: 
Butanoate_metabolism_ko00650<-pathwas_with_KOs[grep("*ko00650*",names(pathwas_with_KOs))][[1]]
Butanoate_metabolism_ko00650_export<-as.data.frame(unique(Butanoate_metabolism_ko00650$KOs), stringsAsFactors=F)
names_genome<-c()
KEGG_genomes_names<-c()
for(i in 1:dim(Butanoate_metabolism_ko00650_export)[1]){
  tmp<-Butanoate_metabolism_ko00650$genome[Butanoate_metabolism_ko00650$KOs == Butanoate_metabolism_ko00650_export[i,1]]
  #adding the name of that genome
  tmp_KEGG<-paste(tmp,collapse = "_")
  tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
  names_genome<-c(names_genome,paste(tmp,collapse = "_"))
  KEGG_genomes_names<-c(KEGG_genomes_names,tmp_KEGG)
}
Butanoate_metabolism_ko00650_export$genomes<-names_genome
Butanoate_metabolism_ko00650_export$genomes<-as.factor(Butanoate_metabolism_ko00650_export$genomes)
levels(Butanoate_metabolism_ko00650_export$genomes) 
Butanoate_metabolism_ko00650_export$genomes_save<-Butanoate_metabolism_ko00650_export$genomes

#we'll produce two maps, with ASD and without ASD
tmp_list<-gsub("ESV","",levels(Butanoate_metabolism_ko00650_export$genomes))
tmp_list<-strsplit(tmp_list,"_")
levels(Butanoate_metabolism_ko00650_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]<-rep("#7733ff",length(levels(Butanoate_metabolism_ko00650_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]))
#Now the ESV with 1-> 10
tmp_list<-gsub("ESV","",levels(Butanoate_metabolism_ko00650_export$genomes))
tmp_list<-as.numeric(tmp_list) #it's ok if there is a warning 
#We need to find exactly the colors for each ESV 
Butanoate_metabolism_ko00650_export$genomes_colors<-Butanoate_metabolism_ko00650_export$genomes
levels(Butanoate_metabolism_ko00650_export$genomes) #"ESV1"    "#7733ff" "ESV15"   "ESV16"   "ESV6"  let's put the colors
levels(Butanoate_metabolism_ko00650_export$genomes)<-c("#f50808","#7733ff","#33c1ff","#183ac6","#f54908")
levels(Butanoate_metabolism_ko00650_export$genomes_colors)<-c("red","purple","light_blue","dark_blue","orange")
colnames(Butanoate_metabolism_ko00650_export)<-c()
write.table(Butanoate_metabolism_ko00650_export, file = "pathway_enrichment/Butanoate_metabolism_ko00650/Butanoate_metabolism_ko00650_export.tsv", quote = FALSE,sep = "\t", row.names =F )

#all the pawtahys are there except lysine degradation which plug into beutanoate patwhays 
#I'm just grabing it ko00310  

#####################################################################################
#####################################################################################
#Now let's look at the pyruvate patwhay
#list wiht all the KOs: list_KO_per_genome
#now th elist of KOs in the lysine pathway: 
lysne_deg_pathway<-c("K00290","K14157","K00293","K14085","K00825","K00164","K00658","K00252","K01692","K01825","K01782",
"K07515","K07514","K07511","K00022","K00626","K18201","K18202","K01582","K14268","K00135","K03897",
"K03896","K03894","K03895","K01843","K01844","K18011","K18012","K18013","K18014","K01034","K01035",
"K21672","K00824","K13609","K19743","K00306","K18854","K06101","K11427","K11420","K09186","K09187",
"K09188","K14959","K09189","K15588","K11424","K11425","K11422","K11423","K19199","K11431","K11428",
"K11421","K18494","K11433","K11419","K11429","K17451","K11430","K11432","K20795","K20796","K18804",
"K18826","K00474","K00128","K00149","K00471","K00473","K13645","K13646","K13647","K11703")

pathwas_with_KOs_lysine<-c()
  for (j in 1:length(list_KO_per_genome)){
    tmp_g<-lysne_deg_pathway[lysne_deg_pathway %in% list_KO_per_genome[[j]]]
    tmp_g<-as.data.frame(tmp_g, stringsAsFactors=F)
    tmp_g$genome_name<-rep(names(list_KO_per_genome)[j],dim(tmp_g)[1])
    colnames(tmp_g)<-c("KOs","genome")
    pathwas_with_KOs_lysine<-rbind(pathwas_with_KOs_lysine,tmp_g)
  }
#and now find for each of them which one is the one in the pathway
names_genome_lys<-c()
KEGG_genomes_names_lys<-c()
for(i in 1:dim(pathwas_with_KOs_lysine)[1]){
  tmp<-Butanoate_metabolism_ko00650$genome[pathwas_with_KOs_lysine$KOs == pathwas_with_KOs_lysine[i,1]]
  #adding the name of that genome
  tmp_KEGG<-paste(tmp,collapse = "_")
  tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
  names_genome_lys<-c(names_genome_lys,paste(tmp,collapse = "_"))
  KEGG_genomes_names_lys<-c(KEGG_genomes_names_lys,tmp_KEGG)
} #all the enzyme are present in both biomarkers types 

####I could dig deeper here but we only have a few genomes sequenced so I'm stopping here 
##########################################################################################################################################
##########################################################################################################################################
#Glycolysis / Gluconeogenesisko00010
Butanoate_metabolism_ko00010<-pathwas_with_KOs[grep("*ko00010*",names(pathwas_with_KOs))][[1]]
Butanoate_metabolism_ko00010_export<-as.data.frame(unique(Butanoate_metabolism_ko00010$KOs), stringsAsFactors=F)
names_genome<-c()
KEGG_genomes_names<-c()
for(i in 1:dim(Butanoate_metabolism_ko00010_export)[1]){
  tmp<-Butanoate_metabolism_ko00010$genome[Butanoate_metabolism_ko00010$KOs == Butanoate_metabolism_ko00010_export[i,1]]
  #adding the name of that genome
  tmp_KEGG<-paste(tmp,collapse = "_")
  tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
  names_genome<-c(names_genome,paste(tmp,collapse = "_"))
  KEGG_genomes_names<-c(KEGG_genomes_names,tmp_KEGG)
}
Butanoate_metabolism_ko00010_export$genomes<-names_genome
Butanoate_metabolism_ko00010_export$genomes<-as.factor(Butanoate_metabolism_ko00010_export$genomes)
levels(Butanoate_metabolism_ko00010_export$genomes) 
Butanoate_metabolism_ko00010_export$genomes_save<-Butanoate_metabolism_ko00010_export$genomes

#we'll produce two maps, with ASD and without ASD
tmp_list<-gsub("ESV","",levels(Butanoate_metabolism_ko00010_export$genomes))
tmp_list<-strsplit(tmp_list,"_")
levels(Butanoate_metabolism_ko00010_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]<-rep("#b1b6bc",length(levels(Butanoate_metabolism_ko00010_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]))
#Now the ESV with 1-> 10
tmp_list<-gsub("ESV","",levels(Butanoate_metabolism_ko00010_export$genomes))
tmp_list<-as.numeric(tmp_list) #it's ok if there is a warning 
levels(Butanoate_metabolism_ko00010_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]<-rep("#f93a3a",length(levels(Butanoate_metabolism_ko00010_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]))
levels(Butanoate_metabolism_ko00010_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]<-rep("#529af0",length(levels(Butanoate_metabolism_ko00010_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]))

colnames(Butanoate_metabolism_ko00010_export)<-c()
write.table(Butanoate_metabolism_ko00010_export, file = "pathway_enrichment/Glycolysis_metabolism_ko00010/Butanoate_metabolism_ko00010_export.tsv", quote = FALSE,sep = "\t", row.names =F )
##########################################################################################################################################
#Glycolysis / Gluconeogenesisko00010
pyruvate_metabolism_ko00620<-pathwas_with_KOs[grep("*ko00620*",names(pathwas_with_KOs))][[1]]
pyruvate_metabolism_ko00620_export<-as.data.frame(unique(pyruvate_metabolism_ko00620$KOs), stringsAsFactors=F)
names_genome<-c()
KEGG_genomes_names<-c()
for(i in 1:dim(pyruvate_metabolism_ko00620_export)[1]){
  tmp<-pyruvate_metabolism_ko00620$genome[pyruvate_metabolism_ko00620$KOs == pyruvate_metabolism_ko00620_export[i,1]]
  #adding the name of that genome
  tmp_KEGG<-paste(tmp,collapse = "_")
  tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
  names_genome<-c(names_genome,paste(tmp,collapse = "_"))
  KEGG_genomes_names<-c(KEGG_genomes_names,tmp_KEGG)
}
pyruvate_metabolism_ko00620_export$genomes<-names_genome
pyruvate_metabolism_ko00620_export$genomes<-as.factor(pyruvate_metabolism_ko00620_export$genomes)
levels(pyruvate_metabolism_ko00620_export$genomes) 
pyruvate_metabolism_ko00620_export$genomes_save<-pyruvate_metabolism_ko00620_export$genomes

#we'll produce two maps, with ASD and without ASD
tmp_list<-gsub("ESV","",levels(pyruvate_metabolism_ko00620_export$genomes))
tmp_list<-strsplit(tmp_list,"_")
levels(pyruvate_metabolism_ko00620_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]<-rep("#b1b6bc",length(levels(pyruvate_metabolism_ko00620_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]))
#Now the ESV with 1-> 10
tmp_list<-gsub("ESV","",levels(pyruvate_metabolism_ko00620_export$genomes))
tmp_list<-as.numeric(tmp_list) #it's ok if there is a warning 
levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]<-rep("#f93a3a",length(levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]))
levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]<-rep("#529af0",length(levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]))
colnames(pyruvate_metabolism_ko00620_export)<-c()
write.table(pyruvate_metabolism_ko00620_export, file = "pathway_enrichment/pyruvate_metabolism_ko00620/pyruvate_metabolism_ko00620_export.tsv", quote = FALSE,sep = "\t", row.names =F )

##########################################################################################################################################
#Tired of doing this, here is a function for it 
for_KEGG<-function(path_ok_grab){
  pyruvate_metabolism_ko00620<-pathwas_with_KOs[grep(path_ok_grab,names(pathwas_with_KOs))][[1]]
  pyruvate_metabolism_ko00620_export<-as.data.frame(unique(pyruvate_metabolism_ko00620$KOs), stringsAsFactors=F)
  names_genome<-c()
  KEGG_genomes_names<-c()
  for(i in 1:dim(pyruvate_metabolism_ko00620_export)[1]){
    tmp<-pyruvate_metabolism_ko00620$genome[pyruvate_metabolism_ko00620$KOs == pyruvate_metabolism_ko00620_export[i,1]]
    #adding the name of that genome
    tmp_KEGG<-paste(tmp,collapse = "_")
    tmp<-genomes_imp$ESV_names[genomes_imp$closest_KEGG_genome %in% tmp]
    names_genome<-c(names_genome,paste(tmp,collapse = "_"))
    KEGG_genomes_names<-c(KEGG_genomes_names,tmp_KEGG)
  }
  pyruvate_metabolism_ko00620_export$genomes<-names_genome
  pyruvate_metabolism_ko00620_export$genomes<-as.factor(pyruvate_metabolism_ko00620_export$genomes)
  levels(pyruvate_metabolism_ko00620_export$genomes) 
  pyruvate_metabolism_ko00620_export$genomes_save<-pyruvate_metabolism_ko00620_export$genomes

#we'll produce two maps, with ASD and without ASD
  tmp_list<-gsub("ESV","",levels(pyruvate_metabolism_ko00620_export$genomes))
  tmp_list<-strsplit(tmp_list,"_")
  levels(pyruvate_metabolism_ko00620_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]<-rep("#b1b6bc",length(levels(pyruvate_metabolism_ko00620_export$genomes)[sapply(tmp_list, function(x) {length(x)>1})]))
#Now the ESV with 1-> 10
  tmp_list<-gsub("ESV","",levels(pyruvate_metabolism_ko00620_export$genomes))
  tmp_list<-as.numeric(tmp_list) #it's ok if there is a warning 
  levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]<-rep("#f93a3a",length(levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>0 & tmp_list<11 & !is.na(tmp_list)]))
  levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]<-rep("#529af0",length(levels(pyruvate_metabolism_ko00620_export$genomes)[tmp_list>10 & tmp_list<22 & !is.na(tmp_list)]))
  colnames(pyruvate_metabolism_ko00620_export)<-c()
  return(pyruvate_metabolism_ko00620_export)
}

#testing the function
propanoate_metabolism<-for_KEGG("*ko00640*")#little edit by hand for color for two ESV same cohort
write.table(propanoate_metabolism, file = "pathway_enrichment/propanoate_metabolism_ko000640/propanoate_metabolism.tsv", quote = FALSE,sep = "\t", row.names =F )

#ok now Sulure metabolilsm 
Sulure_metabolism<-for_KEGG("*920*")
write.table(Sulure_metabolism, file = "pathway_enrichment/Sulure_metabolism_ko000920/Sulure_metabolism.tsv", quote = FALSE,sep = "\t", row.names =F )

deg_aromiatic_ko01220<-for_KEGG("*ko01220*")
write.table(deg_aromiatic_ko01220, file="pathway_enrichment/metabolism_aromatic/deg_aromatic_ko01220.tsv",quote = FALSE,sep = "\t", row.names =F )

flagellar_stuffs<-for_KEGG("*02040*")
write.table(flagellar_stuffs, file="pathway_enrichment/metabolism_flagellar_stuffs_02040/flagellar_stuffs.tsv",quote = FALSE,sep = "\t", row.names =F )

#Aminoacyl-tRNA biosynthesis ko00970
tRNA_stuffs<-for_KEGG("*ko00970*")
write.table(tRNA_stuffs, file="pathway_enrichment/tRNA_stuffs_00970/flagellar_stuffs.tsv",quote = FALSE,sep = "\t", row.names =F )

#Phosphotransferase = transporter ko02060
Phospo_stuffs<-for_KEGG("*ko02060*")
write.table(Phospo_stuffs, file="pathway_enrichment/Phospo_ko02060/Phospo__stuffs.tsv",quote = FALSE,sep = "\t", row.names =F )

#and now for Staphylococcus aureus infection ko05150
Staff_aureus<-for_KEGG("path:ko05150")
write.table(Phospo_stuffs, file="pathway_enrichment/tsaff_aures.tsv",quote = FALSE,sep = "\t", row.names =F )
#methane 
methahne_pathway<-for_KEGG("path:ko00680")
write.table(methahne_pathway, file="pathway_enrichment/methanemethane.tsv",quote = FALSE,sep = "\t", row.names =F )

#metabolites in diverse environment ko01120
diverse_envt_metabolism<-for_KEGG("*ko01120*")#little edit by hand for color for two ESV same cohort
write.table(diverse_envt_metabolism, file = "pathway_enrichment/diverse_envt_metabolism_ko01120/diverse_envt_metabolism.tsv", quote = FALSE,sep = "\t", row.names =F )

##########################################################################################################################################################
##########################################################################################################################################################

##And now let's try to find common metabolites beetween specific pathways
library(bitops)
library(RCurl)
#patwhays: 

path_inte<-readRDS("pathway_enrichment/relevant_GSEA_110.rds")
pathways_imp<-rownames(path_inte)
pathways_imp<-gsub("path:","",pathways_imp)

#basic API: #http://rest.kegg.jp/get/ko00362
grabbing_compounds_in_kegg<-function(pathways_imp){
  path_ESV_choice<-list()
  for (i in 1:length(pathways_imp)){
    cat(i, "\t",pathways_imp[i],"\n")
    path <-getURL (paste("http://rest.kegg.jp/get/",pathways_imp[i],sep=""))
    path<-strsplit(path,"\n")
    path<-path[[1]]
    if (length(grep("*COMPOUND*",path)) > 0 &  length(grep("*REFERENCE*",path)) > 0){
      a<-grep("*COMPOUND*",path)
      b<-min(grep("*REFERENCE*",path))
      path<-path[a:(b-1)]
      path<-gsub("COMPOUND","",path)
      path<-gsub("            ","",path)
      path<-gsub("    ","",path)
      path_ESV_choice[[i]]<-path}
    else if (length(grep("*COMPOUND*",path)) > 0 & length(grep("*REFERENCE*",path)) == 0){
      a<-grep("*COMPOUND*",path)
      b<-length(path)
      path<-path[a:(b-1)]
      path<-gsub("COMPOUND","",path)
      path<-gsub("            ","",path)
      path<-gsub("    ","",path)
      path_ESV_choice[[i]]<-path}
      else
      {path_ESV_choice[[i]]<-path[1]}}
  return(path_ESV_choice)
}


compounds_ok<-grabbing_compounds_in_kegg(pathways_imp)
names(compounds_ok)<-pathways_imp

#now looking for overlap between butanoate and benzoate 
test1<-intersect(compounds_ok[names(compounds_ok) == "ko00362"][[1]],compounds_ok[names(compounds_ok) == "ko00650"][[1]])
test1<-strsplit(test1, "  ")
test1<-sapply(test1, function(x) {x[1]})
cat(test1,sep="\n")
#now looking at butanoate and ko00010= glycolysis 
test1<-intersect(compounds_ok[names(compounds_ok) == "ko00640"][[1]],compounds_ok[names(compounds_ok) == "ko00650"][[1]])
test1<-strsplit(test1, "  ")
test1<-sapply(test1, function(x) {x[1]})
cat(test1,sep="\n")

write.table(test1, file="pathway_enrichment/intersect_benzoate_butanoate.tsv",quote = FALSE,sep = "\t", row.names =F )

          