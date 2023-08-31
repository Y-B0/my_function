name_padding<-function(x,c,f=2,padding=NA){
  ##Automatically populate nonexistent row or column names
  ##x is datafram.
  ##c is name vector
  library(dplyr)
  x<-as.data.frame(x)
  if (f==2) {
    x<-x[,intersect(colnames(x),c),drop=T]
    x[,setdiff(c,colnames(x))]<-padding
    x<-x[,c]
  }else if (f==1) {
    x<-x[intersect(rownames(x),c),,drop=T]
    x[setdiff(c,rownames(x)),]<-padding
    x<-x[,c]
  }
}