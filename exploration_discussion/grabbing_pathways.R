#grabbing the patwhays unique to the genome considered
library(bitops)
library(RCurl)
setwd("Lab/Microbiome_ASD_16S/")

####function to grab pathways form KOs
#Now grabing it in KEGG
grabbing_pathways_in_kegg<-function(ESV_choice_only_KOs){
  path_ESV_choice<-matrix(ncol=3,nrow=length(ESV_choice_only_KOs))
  for (i in 1:length(ESV_choice_only_KOs)){
    cat(i, "\t",ESV_choice_only_KOs[i],"\n")
    path <-getURL (paste("http://rest.kegg.jp/link/path/",ESV_choice_only_KOs[i],sep=""))
    path<-strsplit(path,"\n")
    path<-strsplit(path[[1]],"\t")
    path<-sapply(path, function(x) {x[2]})
    path<-grep("path:ko",path,value=T )
    path<-gsub("path:","",path)
    path_ESV_choice[i,1]<-ESV_choice_only_KOs[i]
    path_ESV_choice[i,2]<-paste(path,collapse = "_")
    path_tmp_s<-c()
    for (j in 1:length(path)){
      path_tmp <-getURL (paste("http://rest.kegg.jp/list/",path[j],sep=""))
      path_tmp<-strsplit(path_tmp,"\t")
      path_tmp_s<-c(path_tmp_s,path_tmp[[1]][2])
    }
    path_tmp_s<-paste(path_tmp_s,collapse = " _ ")
    path_ESV_choice[i,3]<-path_tmp_s
  }
  return(path_ESV_choice)
}

#for ESV8
ESV8_spec_metabo<-read.csv("exploration_discussion/ESV8/ESV8_specificity.csv", stringsAsFactors = F)
ESV8_only_KOs_pres<-ESV8_spec_metabo$KO[ESV8_spec_metabo$PRES_ABS=="pres"]
ESV8_only_KOs_abs<-ESV8_spec_metabo$KO[ESV8_spec_metabo$PRES_ABS=="abs"]

#for ESV6
ESV6_spec_metabo<-read.csv("exploration_discussion/ESV6/ESV6_specificity.csv", stringsAsFactors = F)
ESV6_only_KOs_pres<-ESV6_spec_metabo$KO[ESV6_spec_metabo$PRES_ABS=="pres"]
ESV6_only_KOs_abs<-ESV6_spec_metabo$KO[ESV6_spec_metabo$PRES_ABS=="pres"]


#then apply this function to whatever ESV and/or KOs list you want to look at
path_ESV8_abs<-grabbing_pathways_in_kegg(ESV8_only_KOs_abs)
path_ESV6_abs<-grabbing_pathways_in_kegg(ESV6_only_KOs_abs)


write.csv(path_ESV8, file="exploration_discussion/ESV8/path_ESV8.csv")
write.csv(path_ESV8_abs, file="exploration_discussion/ESV8/path_ESV8_asb.csv")
write.csv(path_ESV6_abs, file="exploration_discussion/ESV6/path_ESV6_asb.csv")

#Looking at how many pathways were detected and how many KOs per pathways 
path_final<-strsplit(path_ESV6_abs[,3], split=" _ ")
list_paths<-unique(unlist(path_final))
list_paths<-table(unlist(path_final))
sort(list_paths)

#for KEGG to look up the pathways, I want to vizualise where this is on KEGG, I just need ot feed the KOs in the GUI interface
path_ESV8[grep("*Biosynthesis of amino acids*",path_ESV8[,3]),1]
path_ESV8[grep("*Biosynthesis of antibiotics*",path_ESV8[,3]),1]

path_ESV8_abs[grep("*Carbohydrate digestion and absorption*",path_ESV8[,3]),1]
path_ESV6_abs[grep("*Biosynthesis of amino acids*",path_ESV6[,3]),1]



