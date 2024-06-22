tpm_deom <- function(expr,na.xy=T){
  human.hg38.protein_coding <- read.delim("E:/project/database/human.hg38.protein_coding.position", header=FALSE)
  if (is.xy) {
    human.hg38.protein_coding<-human.hg38.protein_coding[!human.hg38.protein_coding$Chr%in%c("chrX","chrY"),]
  }
  human.hg38.protein_coding<-human.hg38.protein_coding[!duplicated(human.hg38.protein_coding$Symbol),]
  rownames(human.hg38.protein_coding)<-human.hg38.protein_coding$Symbol
  human.hg38.protein_coding<-human.hg38.protein_coding[intersect(rownames(expr),rownames(human.hg38.protein_coding)),]
  expr<-expr[intersect(rownames(expr),rownames(human.hg38.protein_coding)),]
  kb<-human.hg38.protein_coding$Length/1000
  
  rpk <- expr / kb
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  return(tpm)
}