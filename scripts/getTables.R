library(phyloseq)

getMapping <- function(mp){
  mapping <- read.delim(mp)
  colnames(mapping)[1] = "SampleID"
  mapping$SampleID <- as.character(mapping$SampleID)
  rownames(mapping)<-mapping$SampleID
  mapping$classifier[is.na(mapping$classifier)] <-    4
  mapping <- mapping[order(mapping$SampleID),]
  mapping$Gender <- gsub(" ", "", as.character(mapping$Gender))
  mapping$Treatment <- gsub("c", "C", mapping$Treatment)
  mapping$Treatment <- gsub(" ", "", mapping$Treatment)
  return(mapping)
}

CSS_norm<-function(ps){
  library(metagenomeSeq)
  ps.metaG<-phyloseq_to_metagenomeSeq(ps)
  p_stat = cumNormStatFast(ps.metaG)
  ps.metaG = cumNorm(ps.metaG, p = p_stat)
  ps.metaG.norm <- MRcounts(ps.metaG, norm = T)
  ps_CSS<-phyloseq(otu_table(ps.metaG.norm, taxa_are_rows = T), sample_data(ps),tax_table(ps))
  return(ps_CSS)
}



#Generating a function to have a complete phyloseq object
#argument1:loc.mapping, where the mapping file is
#argument2:loc.biom, which biom file do you want to use 
#argument 3: loc.tree, which tree you want to use
#argument 4: taxa.table, where is the taxa table
#If there is a problem making the phyloseq object, check that all files (biom, mapping, tree, and taxa) are referring to the same step in processing.
#Ie. If you use a taxa table from dada2 but a tree from greengenes, you will see an error. Instead, use the taxa table from qiime with the tree from greengenes

#Suggested usage:
#loc.biom = "ASD_microbiome/data/otu_table_qc_age.biom"
#loc.mapping =  "ASD_microbiome/mapping/mapping_metadata_jan.txt"
#loc.tree = "ASD_microbiome/dada2/tree_unrooted.tre"
#loc.taxa = "ASD_microbiome/dada2/taxa.txt"
#OR
#loc.biom = "ASD_microbiome/data/qiime_otuTables_13.5/otu_table_qc_age.biom"
#loc.mapping =  "ASD_microbiome/mapping/mapping_metadata_jan.txt"
#loc.tree = "ASD_microbiome/data/qiime_otuTables_13.5/gg_13_5_otus_99_annotated.tree"
#loc.taxa = "ASD_microbiome/data/qiime_otuTables_13.5/taxa_onlyThosePresent.txt"


CreatePhyloSeq <- function(loc.mapping,loc.biom,loc.tree, loc.taxa){
  print(loc.biom)
  biom <-phyloseq::import_biom(loc.biom)
  biom <- otu_table(biom, taxa_are_rows = T)
  mapping <- getMapping(loc.mapping)
  
  #Take out all samples from mapping that are not in biom, and vice versa
  print(mapping$SampleID[!(mapping$SampleID %in% colnames(biom))])
  print(colnames(biom) [!(colnames(biom) %in% mapping$SampleID)])
  mapping <- mapping[!(mapping$SampleID %in% setdiff(mapping$SampleID, colnames(biom))),]
  biom <- biom[, !(colnames(biom) %in% setdiff(colnames(biom), mapping$SampleID))]
 
  
  tree<-read_tree_greengenes(loc.tree) # note if it is a green genes tree use function read_tree_greengenes()
 
    
  taxa<-read.csv(loc.taxa, row.names=1)
  
  ps <- phyloseq(otu_table(biom, taxa_are_rows=T), 
                 sample_data(mapping), 
                 tax_table(as.matrix(taxa)), 
                 phy_tree(tree))
  return(ps)
}


CreatePhyloSeq2 <- function(mapping,loc.biom,loc.tree, loc.taxa){
  print(loc.biom)
  biom <-phyloseq::import_biom(loc.biom)
  biom <- otu_table(biom, taxa_are_rows = T)
  
  #Take out all samples from mapping that are not in biom, and vice versa
  mapping <- mapping[!(mapping$SampleID %in% setdiff(mapping$SampleID, colnames(biom))),]
  biom <- biom[, !(colnames(biom) %in% setdiff(colnames(biom), mapping$SampleID))]
  
  
 # tree<-read_tree_greengenes(loc.tree) # note if it is a green genes tree use function read_tree_greengenes()
  
  
  taxa<-read.csv(loc.taxa, row.names=1)
  
  ps <- phyloseq(otu_table(biom, taxa_are_rows=T), 
                 sample_data(mapping), 
                 tax_table(as.matrix(taxa)))
                # phy_tree(tree))
  
  return(ps)
}



#Couldn't figure out how to get phyloseq object to use a given phylogenetic level
#Wrote my own function to condense an otu table by phylogenetic level
#Usage:
#otu_table <- GetOtuTable_byPhylogeneticLevel(phyloseq_obj, "Phylogenetic Level")

#Note:
#Once you do this, you cannot put the table back into a phyloseq object because the taxa names will not match with the taxa table or tree taxa names. 
#You can, however, use for any analysis

GetOtuTable_byPhylogeneticLevel <- function(ps, phyLevel){
    df <- as.data.frame(otu_table(ps))
    levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    level <- match(phyLevel, levels)
    taxaNames <- tax_table(ps)[, levels[1:level]]
    concat_taxaNames <- as.character(apply(taxaNames, 1, paste, collapse=' '))
    df$phy <- concat_taxaNames
    t <- aggregate( .~ phy, df, function(x) sum(x))
    rownames(t) <- t[,1]
    t <- t[,-1]
    return(t)
}


splitAutControl <- function(table, mapping){
  #note that this function deletes any samples that don't have a matching sibling. This
  #is because the analyses that depend on this function assume aut/control tables perfectly align sibling pairs
  
  #KEEP ONLY SIBLINGS
  pairs_to_throw <- names(table(mapping$Pair))[table(mapping$Pair) < 2] #misleading, table is a function to count occurences of pair numbers
  table <- table[, !(mapping$Pair %in% pairs_to_throw)]
  mapping <- mapping[!(mapping$Pair %in% pairs_to_throw), ]
  
  table <- table[ , mapping$SampleID != "45.1"] #thired sibling of pair 2
  mapping <- mapping[ mapping$SampleID != "45.1", ]
  table <- table[ , mapping$SampleID != "233"] #third sibling of pair 55 w/ autism
  mapping <- mapping[ mapping$SampleID != "233", ]
  nsamples <- nrow(mapping)
  
  aut <- table[ , mapping$Treatment == "Aut"]
  control <- table[ , mapping$Treatment == "Control"]
  autMap <- mapping[mapping$Treatment == "Aut",]
  cMap <- mapping[mapping$Treatment == "Control",]
  
  aut <- aut[, order(autMap$Pair)]
  control <- control[,order(cMap$Pair)]
  autMap <- autMap[order(autMap$Pair),]
  cMap <- cMap[order(cMap$Pair),]
  
  ####pairing###
  aut <- aut[ ,autMap$Pair %in% cMap$Pair]
  autMap <- autMap[autMap$Pair %in% cMap$Pair,]
  control <- control[, cMap$Pair %in% autMap$Pair]
  cMap <- cMap[cMap$Pair %in% autMap$Pair,]
  
  return(list(aut, control, autMap,cMap))
}



splitAutControl3_shuff <- function(table, mapping){
  
  #KEEP ONLY SIBLINGS
  pairs_to_throw <- names(table(mapping$Pair))[table(mapping$Pair) < 2] #misleading, table is a function to count occurences of pair numbers
  table <- table[, !(mapping$Pair %in% pairs_to_throw)]
  mapping <- mapping[!(mapping$Pair %in% pairs_to_throw), ]
 
  
  table <- table[ , mapping$SampleID != "45.1"] #thired sibling of pair 2
  mapping <- mapping[ mapping$SampleID != "45.1", ]
  table <- table[ , mapping$SampleID != "233"] #third sibling of pair 55 w/ autism
  mapping <- mapping[ mapping$SampleID != "233", ]
  nsamples <- nrow(mapping)
  

  treatmentVector <- c(rep("Aut", nsamples/2) , rep("Control", nsamples/2))
  mapping$Treatment<-sample(treatmentVector, nsamples,replace = F)
  
  aut <- table[ , mapping$Treatment == "Aut"]
  control <- table[ , mapping$Treatment == "Control"]
  autMap <- mapping[mapping$Treatment == "Aut",]
  cMap <- mapping[mapping$Treatment == "Control",]

  
  #If one has one more sample
 # toRemove = sample(1:max(ncol(aut), ncol(control)), 1)
#  if(ncol(aut) > ncol(control)){
#    aut <- aut[,-toRemove] # get ride of random column
#   autMap <- autMap[-toRemove,]
#  }else if (ncol(control) > ncol(aut)){
#    control <- control[,-toRemove]
#    cMap <- cMap[-toRemove, ]
#  }
  
  autMap$Pair <- sample(seq(1:ncol(aut)), ncol(aut), replace = F)
  cMap$Pair <- sample(seq(1:ncol(control)), ncol(control), replace = F)
  
  aut <- aut[,order(autMap$Pair)]
  control <- control[,order(cMap$Pair)]
  autMap <- autMap[order(autMap$Pair), ]
  cMap <- cMap[order(cMap$Pair), ]
  
  return(list(aut, control, autMap,cMap))
}

#Use if table comes from the biom convert util in qiime
getKoTable <- function(loc.koTxt){
  table <- read.delim(loc.koTxt, row.names = 1)
  colnames(table) <- gsub("X", "", colnames(table))
  table <- table[, order(colnames(table))]
  return(as.data.frame(table))
}

#Use if table comes written manually from R
getKoTable2 <- function(loc.koTxt){
  table <- read.csv("~/Lab/ASD_microbiome/data/piphillan/ko_table_normCSS_piph.txt", row.names=1)
  colnames(table) <- gsub("X", "", colnames(table))
  table <- table[, order(colnames(table))]
  return(as.data.frame(table))
}

#And now plots 

make_boxplots<-function(rsv_table, grouping, extraFeatureLabeling, pvals, title = "", xlab = "", ylab = ""){
  # rsv_table should be a sample x feature table
  # grouping should match the samples and describe what condition the samples have
  # extraFeatureLabeling is a vector that matched the features in the rsv_table that describes extra info you want to display, in this case genus
  #pvals <- round(pvals, digits = 2)
  pvals[pvals < .01] = "***"
  df_grouped <- data.frame(cbind(rsv_table), grouping) # Create dataframe with RSV abundances and Group
  
  colnames(df_grouped) <- c(colnames(rsv_table), "Group") # Rename the columns so we can refer to them later
  
  grouped_stacked <- melt(df_grouped, id = "Group") # Put dataframe in form that's easy to use with ggplot
  
  # Include Genus name in dataframe for graph labelling
  match_seq_to_extraInfo <- data.frame(rsv_name =colnames(rsv_table), extraInfo = extraFeatureLabeling) # Create little mapping dataframe for rsv_names to their genuses
  match_seq_to_pval <- data.frame(rsv_name = colnames(rsv_table), pval_adj = pvals) # Create little mapping dataframe for rsv_names to their genuses
  grouped_stacked$extraInfo <- as.character(match_seq_to_extraInfo$extraInfo[match(grouped_stacked$variable, match_seq_to_extraInfo$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  grouped_stacked$pval <- as.character(match_seq_to_pval$pval_adj[match(grouped_stacked$variable, match_seq_to_pval$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  
  # Plot! The function facet_wrap will break graphs up by whatever variables you put after the '~' character. In this case, we want to break it up by RSV name AND genus name
  p <- ggplot(grouped_stacked, aes(x=Group, y = value)) + geom_boxplot(aes(fill = Group)) +
    geom_jitter(aes(x = Group, y = value), position=position_jitter(0.2), cex=1.5, color="gray44") + facet_wrap(~ variable + extraInfo + pval, scale = "free") + labs(title = title, x = xlab, y = ylab) + scale_y_log10() + theme_minimal()
  
  print(p)
}



