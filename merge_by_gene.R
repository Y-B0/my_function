merge_by_gene<-function(x,symbol,gene,CNV){
  ## label adjust intervals depending on whether the adjacent regions contain the same genes
  ## x        dataset
  ## sep      which symbol to split multiple genes in one interval.
  ## gene     list of gene
  c1<-list(str_split(gene,symbol))[[1]]
  c1<-lapply(c1,function(x){x[x!=""]})
  c2<-CNV
  a=1
  while (a<length(c1)) {
    i<-a
    while (length(intersect(c1[[a]],c1[[a+1]]))!=0 & ((c2[[a]]-2)*(c2[[a+1]]-2)>0)) {
      a=a+1
    }
    x[i:a,"intersect"]<-i
    a<-a+1
    if (a==length(c1)) {
      x[a,"intersect"]<-a
    }
  }
  return(x)
}
