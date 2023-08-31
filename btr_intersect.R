btr_intersect<-function(a,b,...){
  ## a,b    datasets, must
  ## f    the vars in bt.intersect 
  intersect<-list()
  for (i in names(a)) {
    intersect[[i]]<-bt.intersect(a = a[[i]],b = b[[i]],...)
  }
  return(intersect)
}