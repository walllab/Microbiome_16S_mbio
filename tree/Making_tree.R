library(ape)
library(phyloseq)
library(msa)
library("phangorn")
library("phyloseq")
setwd("~/Lab/Microbiome_ASD_16S/")
#tree <- read_tree("tree/tree_unrooted.tre")
archaeSeq = "gatcgattagcatgctagtcgcacgggtttaggcccgtggcggaagctcagtaacacgtggccaaactaccctgtggacgagaataccctcgggaaactgaggtcaattctcgatacggctctcatgctggagtgcagcgagcc
ggaaatgttctggcgccacaggatgtggctgcggccgattaggtagacggtgaggtaacggctcaccgtgccaataatcggtacgggtcatgagag"
archaeSeq <- toupper(archaeSeq)

ps <- readRDS('data/ps_noDuplicates.rds')
otu_table <- otu_table(ps)
seqs <- rownames(otu_table)
seqs <- c(seqs, archaeSeq)
  
mult <- msa(seqs, method="ClustalW", type="dna", order="input")
  
#The phangorn package is then used to construct a phylogenetic tree. Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.
  
phang.align <- msa::msaConvert(mult, type="phangorn::phyDat")
  
phang.align <- as.phyDat(mult, type="DNA")
#phang.align <- readRDS("~/Lab/16S/dada2/phang_align.rds")
attr(phang.align, 'names') <- seqs
#saveRDS(phang.align, paste(opt$outdir, "/phang_align.rds", sep=""))

dm <- dist.ml(phang.align)
  
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$tip.label <- seqs
fit = pml(treeNJ, data=phang.align)
  
## negative edges length changed to 0!
  
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
#write.tree(fitGTR$tree, "tree/tree_mle_gtr_127_noArchae.tre")
  
blo_100<-bootstrap.pml(fitGTR, bs = 100, optNni=TRUE, multicore = TRUE)
treeNJ_bootstrapped_100<-plotBS(fitGTR$tree,blo_100)

#rooted_tree <- root(fitGTR$tree, outgroup= archaeSeq, resolve.root=T)


#and then root the tree?
treeNJ_bootstrapped_rooted_100 <- root(treeNJ_bootstrapped_100, outgroup= archaeSeq, resolve.root=T)

write.tree(treeNJ_bootstrapped_rooted_100, "tree/tree_mle_gtr_rootArchae_127_BS100.tre")
saveRDS(blo,"tree/blo_100.RDS")
saveRDS(fitGTR,"tree/fitGTR.RDS")
  
