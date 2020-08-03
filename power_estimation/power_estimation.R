#load the package for power calculation
source("https://bioconductor.org/biocLite.R")
biocLite("HMP")
setwd("~/Lab/Microbiome_ASD_16S/")
set.seed(100)
library("HMP")

############################################################################################################
#########################################################################################################
#Functions used: 

#########Filtering Taxa
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

#########Calculating the error I and error II
#####Formatting the samples from phyloseq files
#number of reads per samples for ASD samples total 

cal_error_type_I_and_II<- function(ps_filt, MCSIM){ #input: phyloseq and nb of simulation
  otu_ps<-otu_table(ps_filt)
  otu_ps<-as.data.frame(otu_ps)
  totu_ps<-t(otu_ps)
  totu_ps_Aut<-totu_ps[sample_data(ps_filt)$Treatment == "Aut",]
  totu_ps_Aut<-as.data.frame(totu_ps_Aut)
  totu_ps_Aut<-as.matrix(totu_ps_Aut)

  totu_ps_Cont<-totu_ps[sample_data(ps_filt)$Treatment == "Control",]
  totu_ps_Cont<-as.data.frame(totu_ps_Cont)
  totu_ps_Cont<-as.matrix(totu_ps_Cont)

  #### Get a list of dirichlet-multinomial parameters for the data using DOM
  fit.totu_ps_Cont<-DM.MoM(totu_ps_Cont)
  fit.totu_ps_Aut<-DM.MoM(totu_ps_Aut)
  group.alpha<-rbind(fit.totu_ps_Aut$gamma, fit.totu_ps_Cont$gamma)

  ###Seeting up the groups
  nb_reads_totu_ps_Aut<-rowSums(totu_ps_Aut)
  nb_reads_totu_ps_Cont<-rowSums(totu_ps_Cont)
  group.Nrs <- list(nb_reads_totu_ps_Aut, nb_reads_totu_ps_Cont)

  #####Size and Power for the Several-Sample DM Parameter Test Comparison
  #Please set Monte-Carlo to be at least 1,000 (see HMP manual)
  type_II<-MC.Xdc.statistics(group.Nrs,  numMC = MCSIM,group.alpha) #typeII error 0.000999001
  type_I<-MC.Xdc.statistics(group.Nrs,  numMC =MCSIM,fit.totu_ps_Aut$gamma, type="hnull") #typeI error 0.000999001 
  error_typeI_and_II<-list()
  error_typeI_and_II[[1]]<-type_I
  error_typeI_and_II[[2]]<-type_II
  names(error_typeI_and_II)<-c("error_typeI","error_typeII")
  return(error_typeI_and_II)
  }

#########################################################################################################
#########################################################################################################
#####load the ps object and filter it (not normlized)
ps <- readRDS("data/ps_noDuplicates.RDS")
ps<-readRDS("data/ps_filt_norm_117.RDS")
#Removed samples not used for the analysis 
ID_not_ASD<-sample_data(ps)$SampleID[sample_data(ps)$classifier > 0 & sample_data(ps)$Treatment == "Aut"]
#and find thier siblings 
pair_sibling_rem<-sample_data(ps)$Pair[sample_data(ps)$SampleID %in% ID_not_ASD]
samples_to_remove_not_ASD<-sample_data(ps)$SampleID[sample_data(ps)$Pair %in% pair_sibling_rem]
sample_to_keep<-rownames(sample_data(ps))[!rownames(sample_data(ps)) %in% samples_to_remove_not_ASD]
ps <-prune_samples(sample_to_keep,ps)
#filter 
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
otu_ps<-otu_table(ps_filt)


#########################################################################################################
#########################################################################################################
#Power estimation at the ESV level 
otu_ps<-as.data.frame(otu_ps)
totu_ps<-t(otu_ps)

#ok now we;re testing this 
C.alpha.multinomial(totu_ps) #ok reject since zero ##pvalue = 0 => rejected => Dirichlet-Multinomial distribution
Barchart.data(totu_ps, title = "Taxa Proportions")


#power see function above 
cal_error_type_I_and_II(ps_filt,1000)
#$error_typeI
#[1] 0.000999001
#$error_typeII
#[1] 1

#Now let's double check simulating with several samples
#First determine parameters for the population
otu_ps<-otu_table(ps_filt)
otu_ps<-as.data.frame(otu_ps)
totu_ps<-t(otu_ps)
totu_ps_Aut<-totu_ps[sample_data(ps_filt)$Treatment == "Aut",]
totu_ps_Aut<-as.data.frame(totu_ps_Aut)
totu_ps_Aut<-as.matrix(totu_ps_Aut)

totu_ps_Cont<-totu_ps[sample_data(ps_filt)$Treatment == "Control",]
totu_ps_Cont<-as.data.frame(totu_ps_Cont)
totu_ps_Cont<-as.matrix(totu_ps_Cont)
fit.totu_ps_Cont<-DM.MoM(totu_ps_Cont)
fit.totu_ps_Aut<-DM.MoM(totu_ps_Aut)

group.alpha<-rbind(fit.totu_ps_Aut$gamma, fit.totu_ps_Cont$gamma) #those are my population parameters

#Let's run a simulation for en incresing size of samples
a<-seq(from = 10, to = 70, by =4) 
#if we simulated, the average reads per sample is about: 24513
power_0.05=c()
for (i in 1:length(a)){
  cat(i, "\t", tmp_pval, "\t")
  Nrs_ASD <- rep(24000, a[i])
  Nrs_control<-rep(24000, a[i])
  group.Nrs <- list(Nrs_control, Nrs_ASD) #group.alpha stays the same since this was determine for the pop
  tmp_pval<-MC.Xdc.statistics(group.Nrs,  numMC = 1000,group.alpha)
  power_0.05<- c(power_0.05,tmp_pval)
}
saveRDS(power_0.05, file="power_0.05.Rda") #here I skeptical the power should change as the number of sample increases

#plotting 
plot(a,power_0.05) #basically above 40 we're good. I have a hard time to believe those results

#########################################################################################################
#########################################################################################################
#Power estimation at the genus, out of curiosity 

#estimation at the genus levels: Maude is messing around with stuffs
genus_ps<-tax_glom(ps_filt, taxrank="Order")
otu_genus_ps<-otu_table(genus_ps)
otu_genus_ps<-as.data.frame(otu_genus_ps)
totu_genus_ps<-t(otu_genus_ps)
C.alpha.multinomial(totu_genus_ps) #ok
Barchart.data(totu_genus_ps, title = "Taxa Proportions")
cal_error_type_I_and_II(genus_ps,10) #same results

#########################################################################################################
#########################################################################################################
