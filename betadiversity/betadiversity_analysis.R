
ps <- readRDS("data/ps_filt_norm_117.RDS")


library(vegan)
library(reshape2)
library(ggplot2)
ps_aut <- subset_samples(ps, Treatment == "Aut")
ps_nt <- subset_samples(ps, Treatment == "Control")

#Alpha diversity Shannon
aut_div <- diversity(otu_table(ps_aut), index = "shannon", MARGIN = 2)
nt_div <- diversity(otu_table(ps_nt), index = "shannon", MARGIN = 2)
condition <- c(rep("Aut", length(aut_div)), rep("NT", length(nt_div)))
df <- data.frame(condition, div = c(aut_div, nt_div))
div_shannon_plot <- ggplot(df, aes(x = condition, y = div, fill = condition)) + geom_boxplot() + geom_jitter(width = 0.09) +
  ylab("Shannon Diversity Metric") + ggtitle("Shannon Diveristy in Autism and NT siblings")



#Alpha diversity PD
library(picante)
otu_mat <- as(otu_table(ps), "matrix")
pd_div <- pd(samp = t(otu_mat), tree = phy_tree(ps))
rownames(pd_div) == sample_data(ps)$SampleID
pd_div$condition <- sample_data(ps)$Treatment
div_pd_plot <- ggplot(pd_div, aes(x = condition, y = PD, fill = condition)) + geom_boxplot() + geom_jitter(width = 0.09) +
  ylab("Phylogenetic Diversity Metric") + ggtitle("Phylogenetic Diveristy in Autism and NT siblings")


#variance in alpha diversity Shannon: Bootstrap
aut_vars <- c()
nt_vars <- c()
for(i in 1:1000){
  aut_resamples <- sample(names(aut_div), size = length(aut_div), replace = T)
  nt_resamples <- sample(names(nt_div), size = length(nt_div), replace = T)
  aut_vars <- c(aut_vars, var(aut_div[aut_resamples]))
  nt_vars <- c(nt_vars, var(nt_div[nt_resamples]))
}
condition <- c(rep("Aut", length(aut_vars)), rep("NT", length(nt_vars)))
df <- data.frame(condition, vars = c(aut_vars, nt_vars))
var_shannon_plot <- ggplot(df, aes(x = condition, y = vars, fill = condition)) + geom_boxplot() + 
ylab("Variance in Shannon diversity metrics")+
  ggtitle(paste("Variance of Shannon over 10000 sims \n pvalue = ", formatC(wilcox.test(aut_vars, nt_vars)$p.value, format = "e", digits = 2)) )


#Variance in phylogenetic diversity PD: Bootstrap
aut_div <- pd_div[sample_data(ps)$Treatment == "Aut", "PD"]
nt_div <- pd_div[sample_data(ps)$Treatment == "Control", "PD"]
aut_vars <- c()
nt_vars <- c()
for(i in 1:1000){
  aut_resamples <- sample(1:length(aut_div), size = length(aut_div), replace = T)
  nt_resamples <- sample(1:length(nt_div), size = length(nt_div), replace = T)
  aut_vars <- c(aut_vars, var(aut_div[aut_resamples]))
  nt_vars <- c(nt_vars, var(nt_div[nt_resamples]))
}
condition <- c(rep("Aut", length(aut_vars)), rep("NT", length(nt_vars)))
df <- data.frame(condition, vars = c(aut_vars, nt_vars))
var_pd_plot <- ggplot(df, aes(x = condition, y = vars, fill = condition)) + geom_boxplot() + 
  ylab("Variance in Shannon diversity metrics")+
  ggtitle(paste("Variance of PD over 10000 sims \n pvalue = ", formatC(wilcox.test(aut_vars, nt_vars)$p.value, format = "e", digits = 2)) )





ggplot(df, aes(x = condition, y = div, fill = condition)) + geom_boxplot() + geom_jitter(width = 0.09) +
  ylab("Shannon Diversity Metric") + ggtitle("Shannon Diveristy in Autism and NT siblings")


#Diveristy and GI Shannon
div <- diversity(otu_table(ps), index = "shannon", MARGIN = 2)
sample_data(ps)$shannon_div <- div
gi_shannon_plot <- ggplot(sample_data(ps), aes(x = Treatment, y = shannon_div)) + geom_boxplot() + geom_jitter(width = 0.09, aes(colour = bowel.qual, size = 5)) +
  ylab("Shannon Diversity Metric") + ggtitle("Shannon-Weiner Diversity in Autism and NT siblings")

div_s <- sample_data(ps)$shannon_div
sample_data(ps)$shannon_div_cat <- case_when(div_s <= mean(div_s) - sd(div_s) ~ 'Low Div',
                                             between(div_s, mean(div_s) - sd(div_s), mean(div_s) + sd(div_s)) ~ "Med Div",
                                             div_s > mean(div_s) + sd(div_s) ~ 'High Div')

tbl <- table(sample_data(ps)$bowel.qual, sample_data(ps)$shannon_div_cat)
tbl <- tbl[1:3, ]
fisher.test(tbl)



#Diveristy and GI PD
pd_div$GI <- sample_data(ps)$bowel.qual
gi_pd_plot <- ggplot(pd_div, aes(x = condition, y = PD)) + geom_boxplot() + geom_jitter(width = 0.09, aes(colour = GI, size = 5)) +
  ylab("Phylogenetic Diversity Metric") + ggtitle("Phylogenetic Diveristy in Autism and NT siblings")



#Diversity and bowel frequency
bowel_freq <- as.numeric(as.character(sample_data(ps)$bowel.freq))
sample_data(ps)$bowel.freq.cat <- case_when(bowel_freq == 0 ~ "Sparse",
                                            bowel_freq == 1 ~ "Typical",
                                            bowel_freq > 1 ~ "Frequent")
tbl <- table(sample_data(ps)$bowel.freq.cat, sample_data(ps)$shannon_div_cat)
fisher.test(tbl)

#Diversity and Phenotype
#chisq
tbl <- table(sample_data(ps)$Treatment, sample_data(ps)$shannon_div_cat)
fisher.test(tbl)
#wilcoxon
wilcox.test(sample_data(subset_samples(ps, Treatment == "Aut"))$shannon_div,
            sample_data(subset_samples(ps, Treatment == "Control"))$shannon_div)


####PLOT###
library(gridExtra)
grid.arrange(div_pd_plot, var_pd_plot, gi_pd_plot, nrow = 1)
grid.arrange(div_shannon_plot, var_shannon_plot, gi_shannon_plot, nrow = 1)


#############################Other explorations
#beta diversity
aut_dists <- vegdist(t(otu_table(ps_aut)), method = "bray")
nt_dists <- vegdist(t(otu_table(ps_nt)), method = "bray")
boxplot(as.vector(aut_dists), as.vector(nt_dists), col = c("red", "blue"), names = c("Aut", "NT"), ylab = "Within Group Bray Curtis Distances")
legend("topright", fill = c("red", "blue"), legend = c("Autism", "Neurotypical"))

#What is the similarity between the groups that are very close together in ASD
aut_dists <- as.matrix(aut_dists)
nt_dists <- as.matrix(nt_dists)
index <- which(aut_dists < 0.4 & aut_dists > 0, arr.ind = T)
index[,1]

rownames(aut_dists)[index[,1]]
colnames(aut_dists)[index[,2]]
#Wait, they're not actually close to eachother, each just has a very close partner


#Phylum pie charts
ps_aut <- subset_samples(ps, Treatment == "Aut")
tss_aut <- apply(otu_table(ps_aut), 2, function(vec) return(vec/sum(vec) ))
ps_tss_aut <- phyloseq(otu_table(tss_aut, taxa_are_rows = T), tax_table = tax_table(ps), phy_tree = phy_tree(ps), sample_data = sample_data(ps))
ps_tss_aut_glom <- tax_glom(ps_tss_aut, taxrank = "Phylum")

slices_aut <- rowMeans(otu_table(ps_tss_aut_glom))
lbls_aut <- tax_table(ps_tss_aut_glom)[taxa_names(ps_tss_aut_glom), "Phylum"]
pie(slices_aut, labels = lbls_aut, main="Pie Chart of Autism Phylums")

ps_nt <- subset_samples(ps, Treatment == "Control")
tss_nt <- apply(otu_table(ps_nt), 2, function(vec) return(vec /sum(vec) ))
ps_tss_nt <- phyloseq(otu_table(tss_nt, taxa_are_rows = T), tax_table = tax_table(ps), phy_tree = phy_tree(ps), sample_data = sample_data(ps))
ps_tss_nt_glom <- tax_glom(ps_tss_nt, taxrank = "Phylum")

slices_nt <- rowMeans(otu_table(ps_tss_nt_glom))
lbls_nt <- tax_table(ps_tss_nt_glom)[taxa_names(ps_tss_nt_glom), "Phylum"]
pie(slices_nt, labels = lbls_nt, main="Pie Chart of NT Phylums")

cbind(as.vector(slices_aut), as.vector(lbls_aut))
cbind(as.vector(slices_nt), as.vector(lbls_nt))




#PD and dietary supplements
#Diveristy and GI PD
pd_div$supp <- sample_data(ps)$OtherSupp
supp_pd_plot <- ggplot(pd_div, aes(x = supp, y = PD)) + geom_boxplot() + geom_jitter(width = 0.09, aes(colour = supp, size = 5)) +
  ylab("Phylogenetic Diversity Metric") + ggtitle("Phylogenetic Diveristy in Autism and NT siblings")


shannon_div <- data.frame(div = diversity(t(otu_table(ps)), index = "shannon"))
shannon_div$supp <- sample_data(ps)$OtherSupp
supp_shannon_plot <- ggplot(shannon_div, aes(x = supp, y = div)) + geom_boxplot() + geom_jitter(width = 0.09, aes(colour = supp, size = 5)) +
  ylab("Phylogenetic Diversity Metric") + ggtitle("Phylogenetic Diveristy in Autism and NT siblings")


