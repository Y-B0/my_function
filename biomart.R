biomart<-function(x,...){
  library(dplyr)
  library(magrittr)
  library(stringr)
  library(biomaRt)
  
  ## Annotation using genom-site, are
  ## rawdata have colname CHR, start_position, end_position
  ## defult database is human,ensembl hg19
  ## annotated colname are chr, start, end, ensembl, hgnc
  arg<-list(...)
  ensembl<-useMart("ensembl",dataset = "hsapiens_gene_ensembl",host = "grch37.ensembl.org")
  data[,c("chr","start","end","ensembl","hgnc")]<-""
  for (i in 1:nrow(data)) {
    err<-try({
      hg_symbols<- as.data.frame(getBM(attributes=c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol"), filters= c("chromosome_name","start","end"), values = list(chromosome_name=data$CHR[i],start=data$start_position[i],end=data$end_position[i]), mart = ensembl))
    },silent=T
    )
    if (class(err)=="data.frame") {
      hg_symbols<-t(data.frame(apply(err, 2, paste, collapse=",")))
      data[i,c("chr","start","end","ensembl","hgnc")]<-unlist(hg_symbols)
    }else{
      data[i,c("chr","start","end","ensembl","hgnc")]<-"err"
    }
  }
  return(data)
}


