#Here we are performing the analysis for the metadata
library(MASS)
library("ggplot2")
library(pROC)

set.seed(100)
ps<-readRDS("data/ps_filt_norm_117.RDS")
metada<-sample_data(ps)
metada<-as.data.frame(metada)
metada<-as.matrix(metada)
metada<-as.data.frame(metada)

#analysis of GI status: in details: 
metada$bowel.freq.num<-as.numeric(levels(metada$bowel.freq))[metada$bowel.freq]
wilcox.test(bowel.freq.num ~ Treatment, data=metada) #no differences in distribution

#but if we change the category 
metada$bowel.freq.num [metada$bowel.freq.num >1] <- "GI_D"
metada$bowel.freq.num [metada$bowel.freq.num == "0"] <- "GI_D"
metada$bowel.freq.num [metada$bowel.freq.num == "1"] <- "typical"
chisq.test(table(metada$Treatment,metada$bowel.freq.num)) #0.1799 not significant

#Quick terrible plot
metada<-sample_data(ps)
metada<-as.data.frame(metada)
metada<-as.matrix(metada)
metada<-as.data.frame(metada)
metada$bowel.freq.num<-as.numeric(levels(metada$bowel.freq))[metada$bowel.freq]
par(mfrow=c(1,2))
barplot(table(metada$bowel.freq.num[metada$Treatment == "Aut"]), 
        col="red", 
        main="Autism Cohort",
        xlab="Bowel frequency",
        ylab="Count")
barplot(table(metada$bowel.freq.num[metada$Treatment == "Control"]), 
        col="blue",
        main="Neurotypical Cohort",
        xlab="Bowel frequency",
        ylab="Count")

#and last test for this 
library(coin)

#but if we change the category 
metada<-sample_data(ps)
metada<-as.data.frame(metada)
metada<-as.matrix(metada)
metada<-as.data.frame(metada)
metada$bowel.freq.num<-as.numeric(levels(metada$bowel.freq))[metada$bowel.freq]
metada$bowel.freq.num [metada$bowel.freq.num >1] <- "Loose"
metada$bowel.freq.num [metada$bowel.freq.num == "0"] <- "Constipated"
metada$bowel.freq.num [metada$bowel.freq.num == "1"] <- "typical"

tbl_tmp<-table(metada$Treatment,metada$bowel.freq.num)
lbl_test(tbl_tmp) #same results
