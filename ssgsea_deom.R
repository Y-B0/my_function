ssgsea_deom<-function(exp,pathway=NULL,species=="Homo sapiens"){
  library(msigdbr)
  library(GSVA)
  ### exp is a expression data with gene symbol as rowname and sample name as colname
  ### pathway is the gene sets list, if pathway is NULL it will default as kegg for human
  ### species are name in msigdbr
  if (pathway==NULL) {
    pathway <- msigdbr(species = species, category = "C2", subcategory = "KEGG")
    pathway <- pathway %>% split(x = .$gene_symbol, f = .$gs_name)
    names(pathway) <- gsub("KEGG_","",names(pathway))%>%gsub("_"," ",.)
  }

  gs.exp <- gsva(as.matrix(exprSet),pathway, method = "gsva")

}
