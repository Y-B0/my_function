ssgsea_deom<-function(exp,pathway=NULL){
  library(msigdbr)
  library(GSVA)
  ### exp is a expression data with gene symbol as rowname and sample name as colname
  ### pathway is the gene sets list, if pathway is NULL it will default as kegg
  if (pathway==NULL) {
    pathway <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
    pathway <- pathway %>% split(x = .$gene_symbol, f = .$gs_name)
  }

  gs.exp <- gsva(as.matrix(exprSet),pathway, method = "gsva")

}