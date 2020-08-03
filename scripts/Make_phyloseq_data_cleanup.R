###Make Phyloseq
#Removed ids 168.1 and 384 in dada2 quality check
#Removed ids 163 and 211 b/c too old
biom <- import_biom("data/otu_table_qc_age.biom")
otu <- otu_table(biom, taxa_are_rows = T)
map <- readRDS("data/metada.rds")
rownames(map) <- map$SampleID
taxa<- read.csv("dada2/taxa.txt", row.names=1)
tree <- read.tree("tree/tree_mle_gtr_rootArchae_127.tre")
ps <- phyloseq(otu, sample_data(map), tax_table(as.matrix(taxa)), tree)
duplicates <- c("178.6", "179.6", "185.6", "188.6", "220.6", "221.6", "385") # choose to keep 385.6 b/c more reads than its duplicate
ps <- subset_samples(ps, !(SampleID %in% duplicates))
veryLowDiversity <- c("173", "367") #Outliers
ps <- subset_samples(ps, !(SampleID %in% veryLowDiversity))

#saveRDS(ps, "data/ps_noDuplicates.RDS")

write.csv(metada,file="data/metada.cvs")
