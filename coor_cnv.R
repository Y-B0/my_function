library(magrittr)
coor_cnv<-function(n,...){
  ## ...     data_list
  ## n    Specifies the name of each data box to be evaluated by the aspect
  data_list<-list(...)
  names(data_list)<-n
  s_data_list<-data_list
  f_data_list<-data_list
  for (i in 1:length(data_list)) {
    s_data_list[[i]]<-split(data_list[[i]],data_list[[i]]$SAMPLE)
    f_data_list[[i]]<-lapply(s_data_list[[i]], function(x){x<-paste(x$hgnc,collapse=","); x<-list(str_split(x,","))[[1]]; x<-lapply(x,function(x){x[x!=""]})[[1]]})
  }
  intersect_data<-f_data_list[[which(lengths(f_data_list)==min(as.data.frame(lapply(f_data_list,length))))]]
  for (i in names(intersect_data)) {
    #intersect_data[[i]]<-intersect(f_data_list[[1]][[i]],f_data_list[[2]][[i]])
    tmp<-lapply(f_data_list,function(x){x[[i]]})
    intersect_data[[i]]<-Reduce(intersect,tmp)
  }
  intersect_data<-intersect_data[lapply(intersect_data,length)>0]
  s_data_list<-lapply(s_data_list, function(x){x<-x[names(x) %in% names(intersect_data)]})
  for (i in names(intersect_data)) {
    for (j in 1:length(s_data_list)) {
      tmp<-s_data_list[[j]][[i]]
      tmp$hgnc<-as.list(str_split(tmp$hgnc,","))
      tmp<-tmp[lapply(tmp$hgnc, function(x){length(intersect(x,intersect_data[[i]]))})!=0,]
      s_data_list[[j]][[i]]<-tmp
    }
    #s1<-s_data_list[[1]][[i]]
    #s2<-s_data_list[[2]][[i]]
    #s1$hgnc<-as.list(str_split(s1$hgnc,","))
    #s2$hgnc<-as.list(str_split(s2$hgnc,","))
    #s1<-s1[lapply(s1$hgnc, function(x){length(intersect(x,intersect_data[[i]]))})!=0,]
    #s2<-s2[lapply(s2$hgnc, function(x){length(intersect(x,intersect_data[[i]]))})!=0,]
    #s_data_list[[1]][[i]]<-s1
    #s_data_list[[2]][[i]]<-s2
  }
  return(s_data_list)
}
