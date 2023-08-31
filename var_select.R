var_select<-function(x,y){
  ## x: The vector that need to be selected
  rm(list=y[!y%in%c(x,"var_select")])
}