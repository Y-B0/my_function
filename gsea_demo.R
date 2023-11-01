##gsea
gsea<-function(symbol,rank,go.ont="ALL",GO=TRUE,species.go="org.Hs.eg.db",KEGG=TRUE,species.kegg="hsa",n=10){
  ### gene symbol vector for the enrichment, must matched with rank
  ### rank vector for the symbol, only matched with symbol, can be logfc or some sort index for symbol, no need to sort
  ### go.not can be "ALL", "BP", "MF" and "CC"
  ###
  
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(Hmisc)
  library(enrichplot)
  library(RColorBrewer)

  gene <- bitr(symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  deg <- data.frame(gene=gene,rank=rank)
  deg <- dplyr::distinct(deg,gene,.keep_all=TRUE)
  deg <- deg[order(deg$rank,decreasing = T),]
  geneList = deg$rank
  names(geneList) <- deg$ENTREZID
  
    if (GO==T) {
    gse.GO <- try({
      gseGO(
        geneList, 
        ont = go.ont,  
        OrgDb = species,
        keyType = "ENTREZID",
        pvalueCutoff = 1,
        pAdjustMethod = "BH"
      )
    })
    p_go<-try({
      gseaplot2(gse.GO,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))
    })
  }

  if (KEGG==T) {
    gse.KEGG <- try({
      gseKEGG(geneList, 
              organism = "hsa", 
              pvalueCutoff = 1,
              pAdjustMethod = "BH")
    })
    p_kegg<-try({gseaplot2(gse.KEGG,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
  }

  try(return(list(gse.GO=gse.GO,gse.KEGG=gse.KEGG,p_go=p_go,p_kegg=p_kegg)))
}


 
