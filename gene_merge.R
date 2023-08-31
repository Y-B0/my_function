gene_merge<-function(x,col.name){
  ## signal the adjust CNVs, before using this, the data need sort by sample and chr.
  ## colname is the vector of start_interval and end_interval
  data<-x
  data[,"sign"]<-0
  for (i in 1:(nrow(data)-1)) {
    if (data[i,col.name[3]]==data[(i+1),col.name[3]] & (data[i,col.name[2]]==data[(i+1),col.name[1]]|(data[i,col.name[2]])==((data[(i+1),col.name[1]])+1))) {
      data[i:(i+1),"sign"]<-1
    }
    if (data[i,col.name[3]]==data[(i+1),col.name[3]] & data[i,col.name[2]] > data[(i+1),col.name[1]]) {
      data[i:(i+1),"sign"]<-2
    }
  }
  return(data)
}
