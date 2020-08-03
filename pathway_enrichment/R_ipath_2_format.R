#In this script we're realizing a better figure in orde to discuss and vizualized the pathways 
#This is juts ot help to generate txt files that I plugged into ipath and KEGG API. 

gs <- readRDS("pathway_enrichment/pathways_gsea.rds") #all the patwhays from KEGG
ps_piph <- readRDS("data/ps_noDuplicates.RDS") #the taxa detected 

table_norm <- readRDS("pathway_enrichment/ko_abund_table_normByPathway.Rds") #normalized pathways
colnames(table_norm) <- sample_data(ps_piph)$SampleID

#Now we have the list of pathways of importance by the GSEA
relevant_GSEA<-readRDS("pathway_enrichment/relevant_GSEA_110.rds")
selected_GSEA<-relevant_GSEA[relevant_GSEA$q.val<0.05,]
enr_ASD<-selected_GSEA[selected_GSEA$stat.mean>0,]
enr_ASD<-rownames(enr_ASD)
enr_Cont<-selected_GSEA[selected_GSEA$stat.mean<0,]
enr_Cont<-rownames(enr_Cont)

#From there we would like each pathways with colors and which one is enriched where 
#Starting with the ASD: in red for them and larger if KO present 
tmp<-gs[names(gs) == enr_ASD[1]]
#tmp<-paste("\n ", tmp[[1]], sep="")
tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
tmp_NOT_in_data<-tmp[[1]][!tmp[[1]] %in% rownames(table_norm)]
tmp_in_data<-paste(tmp_in_data, "#ff0000 W10")
tmp_NOT_in_data<-paste(tmp_NOT_in_data, "#ff0000 W2")
tmp_to_write_1<-c(tmp_in_data, tmp_NOT_in_data)

tmp<-gs[names(gs) == enr_ASD[2]]
#tmp<-paste("\n ", tmp[[1]], sep="")
tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
tmp_NOT_in_data<-tmp[[1]][!tmp[[1]] %in% rownames(table_norm)]
tmp_in_data<-paste(tmp_in_data, "#ff5400 W10")
tmp_NOT_in_data<-paste(tmp_NOT_in_data, "#ff5400 W2")
tmp_to_write_2<-c(tmp_in_data, tmp_NOT_in_data)
write.table(tmp_to_write_ASD, file="pathway_enrichment/ASD.txt", quote = F, row.names = F, col.names = F)


#and now for the siblings 
colors_nlu<-c("#00d8ff","#00bfff","#00a9ff","#0094ff","#0083ff","#006aff","#003bff","#002aff","#0008ff",
              "#2600ff","#3f00ff","#5900ff", "#7700ff","#00f2ff","#00d8ff")
tmp_to_write_Cont_fin=c()
for (i in 1:length(enr_Cont)){
tmp<-gs[names(gs) == enr_Cont[i]]
tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
tmp_NOT_in_data<-tmp[[1]][!tmp[[1]] %in% rownames(table_norm)]
tmp_in_data<-paste(tmp_in_data," ",colors_nlu[i]," W10", sep="")
tmp_NOT_in_data<-paste(tmp_NOT_in_data," ",colors_nlu[i], " W2", sep="")
tmp_to_write_Cont<-c(tmp_in_data, tmp_NOT_in_data)
tmp_to_write_Cont_fin<-c(tmp_to_write_Cont_fin,tmp_to_write_Cont)}

write.table(tmp_to_write_Cont_fin, file="pathway_enrichment/cont.txt", quote = F, row.names = F, col.names = F)

#Ok it sems like there are some conflicts
#do we have KOs present in both? 
Ko_pathways_in_Cont=c()
for (i in 1:length(enr_Cont)){
  tmp<-gs[names(gs) == enr_Cont[i]]
  tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
  Ko_pathways_in_Cont<-c(Ko_pathways_in_Cont,tmp_in_data)
}
Ko_pathways_in_Cont<-unique(Ko_pathways_in_Cont)

Ko_pathways_in_ASD=c()
for (i in 1:length(enr_ASD)){
  tmp<-gs[names(gs) == enr_ASD[i]]
  tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
  Ko_pathways_in_ASD<-c(Ko_pathways_in_ASD,tmp_in_data)
}
Ko_pathways_in_ASD<-unique(Ko_pathways_in_ASD)

KO_in_both<-intersect(Ko_pathways_in_Cont,Ko_pathways_in_ASD)
#"K01885" "K02398" "K02402" "K02403" "K02405" "K02406" "K02556"

#I will extract by hand all the patwhays 
Ko_pathways_in_Cont=c()
for (i in 1:length(enr_Cont)){
  tmp<-gs[names(gs) == enr_Cont[i]]
  tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
  tmp_in_data<-paste(tmp_in_data, "blue")
  write.table(tmp_in_data, file="pathway_enrichment/tmp_in_data.txt",quote = F, row.names = F, col.names = F)

#and the ASD ones 
Ko_pathways_in_ASD=c()
for (i in 1:length(enr_ASD))){
  tmp<-gs[names(gs) == enr_ASD[i]]
  tmp_in_data<-tmp[[1]][tmp[[1]] %in% rownames(table_norm)]
  tmp_in_data<-paste(tmp_in_data, "orange")
  write.table(tmp_in_data, file="pathway_enrichment/tmp_in_data.txt",quote = F, row.names = F, col.names = F)
}
