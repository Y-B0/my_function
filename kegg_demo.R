kegg_demo<-function(genesymbol,ntop=10,plot=T,species=org.Hs.eg.db){
  library(clusterProfiler) 
  library(stringr) 
  library(AnnotationDbi)
  library(org.Hs.eg.db) 
  library(org.Mm.eg.db)
  library(DOSE)
  library(ggplot2) 
  library(ggrepel) 
  library(ggsci)
  
  try({
    id_list <- mapIds(species,genesymbol,"ENTREZID","SYMBOL")
    id_list <- na.omit(id_list)
    
    kegg <- enrichKEGG(id_list, keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = "BH",
                       minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE)
    kegg<-DOSE::setReadable(kegg, OrgDb='org.Hs.eg.db', keyType = 'ENTREZID')
    
    p<-ggplot(kegg[1:ntop,], aes(x=GeneRatio, y=Description,size=Count, color=pvalue)) + geom_point() + 
      scale_colour_gradient(low="green",high="red") + labs(color=expression(padj),size="Gene number", x="GeneRatio",y="Pathway name",title="KEGG Pathway enrichment")
    return(list(kegg=kegg,plot=p))
  })
}
