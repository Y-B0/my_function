##
doimod <- function(rawpath,refpath) {
  library(stringr)
  library(do)
  library(xfun)
  ref<-readLines(refpath)
  ref<-ref[grep("\\@.*\\{.*|doi = ",ref)]
  ref<-trimws(ref)
  i<-2
  while (i<=length(ref)) {
    if (grepl("\\@.*\\{",ref[i])) {
      ref<-append(ref,paste("[@",str_split(ref[i-1],"\\{|,",simplify = T)[2],"]",sep = ""),after = (i-1))
    }
    i<-i+2
  }
  if (i>length(ref)) {
    ref[i]<-paste("[@",str_split(ref[i-1],"\\{|,",simplify = T)[2],"]",sep = "")
  }
  
  ref<-data.frame(index=ref[seq(1,length(ref),2)],doi=ref[seq(2,length(ref),2)])
  ref$index<-Replace(ref$index,pattern = c("\\@.*\\{:[\\@",",:]"),)
  ref$doi<-Replace(ref$doi,pattern = c("doi = \\{:[\\@","\\},:]"))
  key<-paste(ref$doi,ref$index,sep = ":")
  
  apply(ref, 1, function(x){
    gsub_file(rawpath,x[2],x[1],fixed=T)
  })
  
}

