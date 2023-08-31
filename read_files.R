read_files<-function(path,...){
  para<-list(...)
  filename<- list.files(path) 
  filepath<-file.path(path,filename)
  filelist<-list()
  for (i in 1:length(filepath)) {
    try({
      filelist[[i]]<-read.table(filepath[[i]],...)
      names(filelist)[[i]]<-str_split(filepath[[i]],"/",simplify = T)[length(str_split(filepath[[i]],"/",simplify = T))]
      })
  }
  return(filelist)
}
