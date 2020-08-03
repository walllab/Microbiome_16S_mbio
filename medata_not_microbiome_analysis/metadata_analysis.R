#Here we are performing the analysis for the metadata
library(MASS)
library("ggplot2")
library(pROC)

set.seed(100)
ps<-readRDS("data/ps_filt_norm_117.RDS")
metada<-sample_data(ps)
metada<-as.data.frame(metada)


#Double checking all of the 57 questions 
for (i in 1:length(colnames(metada))){
cat(colnames(metada)[i],levels(unlist(metada[,i])),"\n")
}

#ok now looks good. To find the ordinal vs categorical: the ones with more than 2 levels! (beside the numerical ones) 
to_analyze<-which(!colnames(metada)%in%c("SampleID","Treatment","Pair","batch","classifier","Current.Res","Travel","Roommates","sev","severe"))
metada_to_analyze<-metada[,to_analyze]
length_lev<-c()
for (i in 1:length(colnames(metada_to_analyze))){
  length_lev_tmp<-length(levels(unlist(metada_to_analyze[,i])))
  length_lev<-c(length_lev,length_lev_tmp)
}

#Categorical data
categorical_data<-colnames(metada_to_analyze)[length_lev ==2]
ordinal_data<-colnames(metada_to_analyze)[length_lev>2]

####Analysis catergorical data 
p.val_save=c()
for (i in 1:length(categorical_data)){
  tmp<-metada[,categorical_data[i]]
  tmp<-as.list(tmp)
  tbl_tmp<-table(metada$Treatment,tmp[[1]])
  p.val_tmp<-chisq.test(tbl_tmp)$p.value
  names(p.val_tmp)<-categorical_data[i]
  p.val_save<-c(p.val_save,p.val_tmp)
}

p.val_save_adj_chi_square<-p.adjust( p.val_save,method ="fdr")
p.val_save_adj_chi_square<-p.val_save_adj_chi_square[p.val_save_adj_chi_square<0.05]
print(p.val_save_adj_chi_square) 

#ALL
#Gender       OtherSupp other.diet.rest          gluten           dairy 
#0.0233973360    0.0227425626    0.0227425626    0.0003515719    0.0025998473

#Now the ordinal data
library(coin)
p.val_save_ordinal=c()
for (i in 1:length(ordinal_data)){
  tmp<-metada[,ordinal_data[i]]
  tmp<-as.list(tmp)
  tbl_tmp<-table(metada$Treatment,tmp[[1]])
  p.val_tmp<-pvalue(lbl_test(tbl_tmp))
  names(p.val_tmp)<-ordinal_data[i]
  p.val_save_ordinal<-c(p.val_save_ordinal,p.val_tmp)
}
p.val_save_ordinal_adj<-p.adjust( p.val_save_ordinal,method ="fdr")
p.val_save_ordinal_adj<-p.val_save_ordinal_adj[p.val_save_ordinal_adj<0.05]
print(p.val_save_ordinal_adj) 
#Milk..Cheese 
#0.002096976

#and finally just a t-test for the age
#wilcox.test(metada$age[metada$Treatment == "Aut"],metada$age[metada$Treatment == "Control"])
#Not significant 
#I also did a bunch of tests for the stool constitancy and did not find it either 

######################################################################################################################################## 
#reminder functions splitAutControl is in the R file called "scripts.R"
#Ok now we're doing distances ramdon or betweent the different 
# we need ot change all variable with number for the euclidean distances 
metada_to_analyze_num<-as.list(metada_to_analyze)
metada_to_analyze_num<-lapply(metada_to_analyze_num,function(x) as.numeric(x))
metada_to_analyze_num <- matrix(unlist(metada_to_analyze_num), ncol = 47, byrow = FALSE)
colnames(metada_to_analyze_num)<-names(as.list(metada_to_analyze))
rownames(metada_to_analyze_num)<-metada_to_analyze@row.names
tmetada_to_analyze<-t(metada_to_analyze_num)

#ok that we use for the mean 
split_metada<-splitAutControl(tmetada_to_analyze,metada)
#Distance between the pairs form the same family
dis_sib<-c()
for (i in 1:dim(split_metada[[1]])[2]){
tmp<-dist(rbind(split_metada[[1]][,i], split_metada[[2]][,i]))
dis_sib<-c(dis_sib,tmp)
}
print (mean(dis_sib))
dis_sib_final_sabe<-dis_sib
#and now with permutation 
#basically do the same BUT you shuffle first

mean_shuff_eucl=c()
 for (j in 1:999){
metada_shuff<-metada
metada_shuff$Pair<-sample(metada_shuff$Pair)
split_metada<-splitAutControl(tmetada_to_analyze,metada_shuff)
dis_sib<-c()
  for (k in 1:min(dim(split_metada[[1]])[2],dim(split_metada[[2]])[2])){
    tmp<-dist(rbind(split_metada[[1]][,k], split_metada[[2]][,k]))
    dis_sib<-c(dis_sib,tmp)
  }
tmp_mean<-mean(dis_sib)
mean_shuff_eucl<-c(mean_shuff_eucl,tmp_mean)}
mean_shuff_eucl_save<-mean_shuff_eucl
save_data_eucl=list()
save_data_eucl[[1]]<-dis_sib_final_sabe
save_data_eucl[[2]]<-mean_shuff_eucl_save
save(save_data_eucl,file="medata_not_microbiome_analysis/save_data_eucl.Rda")
#and now just calcuate the pval of this 
load("medata_not_microbiome_analysis/save_data_eucl.Rda")
#ok now plot it 
# generate the plot

#calculate the pval : no need it's Inf
aa<-as.data.frame(save_data_eucl[[2]])
p1<-ggplot(aa) + stat_density(aes(save_data_eucl[[2]]),geom="line",position="dodge", colour = "black", size = 0.3 )

pdf("medata_not_microbiome_analysis/Euclidean_distances_sib.pdf",width=8,height=4)
p1 + geom_vline(xintercept = mean(save_data_eucl[[1]]), colour = "yellow", size= 1 ) + 
  labs(y= "Density",x = "Euclidean Distances", title="Euclidean Distances of 999 permutation compared to Euclidean Distances between siblings")
dev.off()

