pre_plot<-function(datasets,col,group="NA",col.name="NA"){
  ##merge mutil dataset to a dataframe to plot
  ##datasets:a list of all dataset
  ##col:a list about datasets colnames, the order must corresponding with datasets
  ##group:a vector that group the final datasets
  ##col.name:new colname, must set, if don't set when colname set to 1:ncol(data)
  datalist<-list()
  if (group!="NA") {
    for (i in 1:length(datasets)) {
      datalist[[i]]<-datasets[[i]][,col[[i]]]
      datalist[[i]][,"group"]<-group[i]
      if (col.name!="NA") {
        colnames(datalist[[i]])<-col.name
      }else{
        colnames(datalist[[i]])<-paste("V",1:ncol(datalist[[i]]),sep = "")
      }
    }
  }else if (length(names(datalist)!=0)) {
    for (i in 1:length(datasets)) {
      datalist[[i]]<-datasets[[i]][,col[[i]]]
      datalist[[i]][,"group"]<-names(datasets)[i]
      if (col.name=="NA") {
        colnames(datalist[[i]])<-paste("V",1:ncol(datalist[[i]]),sep = "")
      }else{
        colnames(datalist[[i]])<-col.name
      }
    }
  }else{
    for (i in 1:length(datasets)) {
      datalist[[i]]<-datasets[[i]][,col[[i]]]
      datalist[[i]][,"group"]<-i
      if (col.name=="NA") {
        colnames(datalist[[i]])<-paste("V",1:ncol(datalist[[i]]),sep = "")
      }else{
        colnames(datalist[[i]])<-col.name
      }
    }
  }
  data<-Reduce(rbind,datalist)
  return(data)
}

