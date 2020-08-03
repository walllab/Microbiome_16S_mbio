#Pathway enrichment GSEA plotting

plotPValues <- function(greater, less, condition=""){
  op <- par(mar = c(2,15,10,0) )
  
  vals_enriched <- sort(p.adjust(greater[,'p.val'], method="fdr"), decreasing = F)
  #barplot( vals_enriched[1:10], horiz = T, las=1, cex.names=0.6, main = paste("Enriched in Autism", "\n", condition, sep= " "), xlab = "P value")
  
  
  print(vals_enriched[ vals_enriched < 0.1])
  vals_enriched_graph <- greater[greater[,'q.val'] < .05, 'q.val']
  vals_enriched_graph <- vals_enriched_graph[!is.na(vals_enriched_graph)]
  
  vals_depleted <- sort(p.adjust(less[,'p.val'], method="fdr"), decreasing = F)
  #barplot(vals_depleted[1:10], horiz = T, las=1, cex.names=0.6, main = paste("Depleted in Autism", "\n", condition, sep= " "), xlab = "P value")
  
  
  print(vals_depleted[ vals_depleted < 0.1])
  vals_depleted_graph <- less[less[,'q.val'] < .05, 'q.val']
  vals_depleted_graph <- vals_depleted_graph[!is.na(vals_depleted_graph)]
  
  
  if((length(vals_enriched_graph) + length(vals_depleted_graph)) > 0){
    df <- data.frame(
      var1 = paste(1:(length(vals_enriched_graph) + length(vals_depleted_graph)), c(names(vals_enriched_graph), names(vals_depleted_graph))),
      var2 = c(-log(vals_enriched_graph), log(vals_depleted_graph)),
      enriched = c(rep('red', length(vals_enriched_graph)), rep("blue", length(vals_depleted_graph)))
    )
    df$var1 <- factor(df$var1, levels = df$var1)
    
    
    par(mar = c(2,20,10,5))
    p <- ggplot(df, aes(x=var1, y= var2)) +
      geom_bar(stat="identity", position="dodge", fill = df$enriched) +
      guides(fill=FALSE)+
      coord_flip()+
      ylab("-log(PValue)") + xlab("KEGG Ortholog")+
      labs(title = condition)
    
    print(p)
  }
  
  print(paste("Aut sample size: ", sum(map$Treatment == "Aut")))
  print(paste("Control sample size: ", sum(map$Treatment == "Control")))
} 
