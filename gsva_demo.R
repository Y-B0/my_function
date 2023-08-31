gsva<-function(gmt.path="~/Documents/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c2.cp.v7.4.symbols.gmt",exprSet,pathway=NULL){
  ##GSVA
  ##gmt.path file :c2.cp.v7.4.symbols.gmt
  ## row means genes and cols means samples
  
  library(GSVA)
  library(GSEABase)
  c2gmt <- getGmt(gmt.path)
  if (!is.null(pathway)) {
    c2gmt <- c2gmt[grep(pathway, names(c2gmt)),]
  }
  gs.exp <- GSVA::gsva(as.matrix(exprSet), c2gmt, kcdf = "Poisson", min.sz = 10, method='gsva',abs.ranking=TRUE)
  return(gs.exp)
}