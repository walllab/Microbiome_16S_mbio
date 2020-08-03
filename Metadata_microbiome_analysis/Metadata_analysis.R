require(vegan)
library(phyloseq)
library(ggplot2)
library(ade4)
library(gridExtra)
library(grid)
library(lattice)
library(ape)

#setwd("/Users/mmdavid/Microbiome_ok/Microbiome_ASD_16S/")
#biom <- "data/otu_table_qc_age.biom"
#loc.mapping = "mapping/mapping_file_May2017_FINAL.txt" ##This mapping file does not have the smaples with low quality
#mapping <- sample_data(getMapping(loc.mapping))
#making sure the number are not considered as value but factors
#mapping$batch<-as.factor(mapping$batch)
#mapping$pair<-as.factor(mapping$pair)

#treeok<- "tree/rootedtree_long_branch.tre"
#loc.taxa <- "dada2/taxa.txt"
#source ("scripts/getTables_christine.R")

#ps<-CreatePhyloSeq2(mapping,biom,treeok,loc.taxa)

#ok here we will directly upload the CSS norm file
#ps.norm.CSS
setwd("~/Lab/Microbiome_ASD_16S/")
set.seed(100)

#################################################Declaring all of our functions#######################
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


#Loading data
load("Metadata_microbiome_analysis/metada.Rda")
ps <- readRDS("data/ps_noDuplicates.RDS")
phy_tree(ps) <- read.tree("tree/tree_mle_gtr_rootArchae_127.tre")
#sample_data(ps)<-metada #already loaded in ps


#Remove the samples which don't have ASD and thier sibling: 
ID_not_ASD<-sample_data(ps)$SampleID[sample_data(ps)$classifier > 0 & sample_data(ps)$Treatment == "Aut"]
#and find thier siblings 
pair_sibling_rem<-sample_data(ps)$Pair[sample_data(ps)$SampleID %in% ID_not_ASD]
samples_to_remove_not_ASD<-sample_data(ps)$SampleID[sample_data(ps)$Pair %in% pair_sibling_rem]
sample_to_keep<-rownames(sample_data(ps))[!rownames(sample_data(ps)) %in% samples_to_remove_not_ASD]
ps<-prune_samples(sample_to_keep,ps)

###Filter and normalize
prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
ps_filt_norm <- deSeqNorm(ps_filt)
saveRDS(ps_filt_norm,file="data/ps_filt_norm_117.RDS")
mapping<-sample_data(ps_filt_norm)
mapping<-as.data.frame(mapping)


#Calculate the distances
#now only keep the factors for each analysis, i.e. remove the samples with no data

pval_factors=c()
nb_samples=c()
for (i in 2:length(mapping)){
  cat (i,"\t")
  tmp_map<-mapping[!is.na(mapping[,i]),]
  tmp_map<-mapping[!mapping[,i] == "",]
  ps.tmp<-ps_filt_norm
  sample_data(ps.tmp) <- tmp_map
  #OTU_tmp<-CreatePhyloSeq2(tmp_map,biom,treeok,loc.taxa)
  tmp_nb_samples<-dim(otu_table(ps.tmp))[2]
  OTU_tables_bray <- phyloseq::distance(ps.tmp, method = "bray")
  df_metadata <- data.frame(sample_data(ps.tmp ))
  form1<-as.formula(paste("OTU_tables_bray",colnames(mapping)[i],sep="~"))
  tmp<-adonis(form1, data = df_metadata, permutations = 9999)
  tmp<-tmp$aov.tab$`Pr(>F)`[1]
  pval_factors<-c(pval_factors,tmp)
  nb_samples<-c(nb_samples,tmp_nb_samples)}

names(pval_factors) <- colnames(mapping)[2:length(mapping)]
pval_factors<-p.adjust(pval_factors, method = "fdr")
save.pval.adonis0.05_braycurtis<-pval_factors[pval_factors<0.05]
save.nbsamples.pva.adonis.05<-nb_samples[pval_factors<0.05]
print(save.pval.adonis0.05_braycurtis)
#print(save.nbsamples.pva.adonis.05)
save(save.pval.adonis0.05_braycurtis,save.nbsamples.pva.adonis.05, file="Metadata_microbiome_analysis/save.pval.adonis.05.Rda")


pval_factors_diper=c()
nb_samples_disper=c()
for (i in 2:length(mapping)){
  cat (i,"\t")
  test_map<-mapping[!is.na(mapping[,i]) & mapping[,i] != "" ,]
  ps.tmp<-ps_filt_norm
  sample_data(ps.tmp) <- test_map
  #OTU_tmp<-CreatePhyloSeq2(test_map,biom,treeok,loc.taxa)
  #OTU_tmp <- subset_samples(OTU_tmp, colnames(mapping)[i] !="")
  df_metadata <- data.frame(sample_data(ps.tmp))
  df_metadata<-df_metadata[df_metadata[,colnames(test_map)[i]] != "",]
  ps.tmp <- subset_samples(ps.tmp, colnames(test_map)[i] !="")
  tmp_nb_samples<-dim(otu_table(ps.tmp))[2]
  OTU_tables_bray <- phyloseq::distance(ps.tmp, method = "bray")
  beta <- betadisper(OTU_tables_bray, df_metadata[,colnames(test_map)[i]])
  tmp<-permutest(beta)
  tmp<-tmp$tab$`Pr(>F)`[1]
  pval_factors_diper<-c(pval_factors_diper,tmp)
  nb_samples_disper<-c(nb_samples_disper,tmp_nb_samples)}

names(pval_factors_diper) <- colnames(mapping)[2:length(mapping)]
pval_factors_diper<-p.adjust(pval_factors_diper, method = "fdr")
pval_factors_diper0.05<-pval_factors_diper[pval_factors_diper>0.05] #those are the one with homogeneous dispersion
nb_sample_disper0.05<-nb_samples_disper[pval_factors_diper>0.05]
print(pval_factors_diper0.05)
#print(nb_sample_disper0.05)

potential_confounding_factor<-intersect(names(pval_factors_diper0.05),names(save.pval.adonis0.05_braycurtis))
print(potential_confounding_factor)
save(pval_factors_diper0.05,nb_sample_disper0.05, file="Metadata_microbiome_analysis/save.pval.adonis.05.Rda")
save(potential_confounding_factor, file="Metadata_microbiome_analysis/potential_confounding_factor.Rda")


####################################################plotting everything#################
#Here are the possible confonding factors: 
#"batch"        "Multivitamin" "Meat.eggs    "Sugar"        "OliveOil" 

load("Metadata_microbiome_analysis/potential_confounding_factor.Rda")
#Needs to remove the NA
conf_map<-as.data.frame(mapping)
conf_map<-mapping[,colnames(mapping) %in% c(potential_confounding_factor)]
conf_map<-conf_map[complete.cases(conf_map), ]
ps.tmp<-ps_filt_norm
sample_data(ps.tmp) <- conf_map

#Looking at those with ordination and constraintes
#ps_not_na <- subset_samples(physeq = ps.norm.CSS,
#  !is.na (Multivitamin) & 
#    !is.na (Frozen.D) &
#    !is.na (Sugar) & 
#    !is.na (Meat.eggs) &
#    !is.na (Water.Day))


ps_pcoa <- ordinate(
  physeq = ps.tmp, 
  method = "CAP", 
  distance = "bray",
  formula = ~ Multivitamin + batch + Meat.eggs + Sugar + OliveOil 
)


# CAP plot, to arrange on a grid after:
#For graphics: the titles:
tite_prep<-c("Panel 1: Bacth sequencing",
             "Panel 2: Is he/she taking a daily multivitamin?",
             "Panel 3: In an average week, how often does your child consume meat/eggs?",
             "Panel 4: How many days in a week does your child consume sugary sweets?",
              "Panel 5: How frequenlty do you cook with olive oil?")

to_plot=list()
for (i in 1:length(potential_confounding_factor)){
  to_plot[[i]] <- plot_ordination(
  physeq = ps.tmp, 
  ordination = ps_pcoa, 
  color = potential_confounding_factor[i], 
  axes = c(1,2),
  title=tite_prep[i]
) + 
  geom_point( size = 2) +
  theme(text = element_text(size =10), plot.title = element_text(size=10))
}

#I want to change the levels just for more visiblity while plotting 
levels(to_plot[[1]]$data$Multivitamin)<-c("yes","No")
levels(to_plot[[2]]$data$batch)<-c("0","1", "2","3","4","5","6")
levels(to_plot[[3]]$data$Meat.eggs)<-c("Never","Rarely(a few times a month)", "Ocasionally(1-2 times a week)","Regularly(3-5/week","Daily")
levels(to_plot[[4]]$data$Sugar)<-c("Never","Rarely(a few times a month)", "Ocasionally(1-2 times a week)","Regularly(3-5/week","Daily")
levels(to_plot[[5]]$data$OliveOil)<-c("Never","Rarely(a few times a month)", "Ocasionally(1-2 times a week)","Regularly(3-5/week","Daily")

#adding the plots with the taxa to be able to 

to_plot[[6]] <-plot_ordination(physeq = ps.tmp, ordination = ps_pcoa, type="taxa", color="Phylum", title ="Panel 7: Taxa") + theme(text = element_text(size =8))
lay <- rbind(c(1,2),c(3,4),c(5,6))

pdf("Metadata_microbiome_analysis/confounding_factors.pdf",width=11,height=11)

#do.call("grid.arrange", c(plot_order_mapped, ncol=2))
#Wait I can do better 
grid.arrange(grobs = to_plot, layout_matrix = lay)
dev.off()
#ok and now all I need is to be able to be able to change the legend in ggplot2

#ok let's try to find the spcies that show some importance in this PCA
taxa.to.select<-vegan::scores(ps_pcoa)$species
#now plot it with no name for visibilty
rownames(taxa.to.select)<-c()
s.arrow(taxa.to.select) #the taxa that influence the most the plots are above 0.25
taxa.to.select.to.rem<-vegan::scores(ps_pcoa)$species[abs(vegan::scores(ps_pcoa)$species[,1])>0.1 | abs(vegan::scores(ps_pcoa)$species[,2])>0.1,]
save(taxa.to.select.to.rem, file="Metadata_microbiome_analysis/taxa.to.select.to.rem.Rda")


save(metada,file="Metadata_microbiome_analysis/metada.Rda")
save_taxa_to_remove<-tax_table(ps)[rownames(tax_table(ps))%in% rownames(taxa.to.select.to.rem),]
save(save_taxa_to_remove,file="Metadata_microbiome_analysis/save_taxa_to_remove.Rda")
write.csv(save_taxa_to_remove,file="Metadata_microbiome_analysis/save_taxa_to_remove.csv",quote=F)


##################################################################################
#IF we want to addd arrow and more stuffs, not used in the paper, free code for whoever wants to have lots of fun
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(ps_pcoa, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.8 * CAP1, 
                 y = 1.8 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
pdf("Metadata_microbiome_analysis/PcoA_Factors_metadata.pdf")
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 3,  
    data = arrowdf, 
    show.legend = FALSE
  )
dev.off()

anova(ps_pcoa) #save this PDF


