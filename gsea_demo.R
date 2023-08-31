##gsea
gsea<-function(deg,fc.name="logFC",go.ont="BP",n=10){
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(Hmisc)
  library(enrichplot)
  library(RColorBrewer)
  gene=bitr(rownames(deg),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  deg<-deg[match(intersect(rownames(deg),gene$SYMBOL),rownames(deg)),]
  gene<-gene[match(intersect(rownames(deg),gene$SYMBOL),gene$SYMBOL),]
  deg$ENTREZID<-gene$ENTREZID
  deg<-deg[order(deg$logFC,decreasing = T),]
  geneList = deg[,fc.name]
  names(geneList) <- deg$ENTREZID
  
  gse.GO <- try({
    gseGO(
      geneList, 
      ont = go.ont,  
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      pvalueCutoff = 1,
      pAdjustMethod = "BH"
    )
  })
  p_go<-try({
    gseaplot2(gse.GO,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))
  })
  #write.csv(gse.GO[[1]],"gse_GO.csv",row.names = F)
  gse.KEGG <- try({
    gseKEGG(geneList, 
            organism = "hsa", 
            pvalueCutoff = 1,
            pAdjustMethod = "BH")
  })
  p_kegg<-try({gseaplot2(gse.KEGG,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
  try(return(list(gse.GO=gse.GO,gse.KEGG=gse.KEGG,p_go=p_go,p_kegg=p_kegg)))
  #write.csv(gse.KEGG[[1]],"gse_KEGG.csv",row.names = F)
}


 