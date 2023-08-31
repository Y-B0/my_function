cor_melt<-function(cormat,pmat) {
  library(reshape2)
  cormat<-melt(cormat)
  pmat<-melt(pmat)
  d<-data.frame(cormat,pmat$value)
  colnames(d)<-c("row","col","cor_value","p_value")
  return(d)
}
