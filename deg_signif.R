deg_signif<-function(dat,signif.name="P.Value",fc.name="logFC",fc.num=0.585,p.num=0.05){
  x<-dat
  x$Sig[(x[,signif.name] > p.num|x[,signif.name]=="NA")|(x[,fc.name] < fc.num)& x[,fc.name] > -fc.num] <- "Not"
  x$Sig[x[,signif.name] <= p.num & x[,fc.name] >= fc.num] <- "Up"
  x$Sig[x[,signif.name] <= p.num & x[,fc.name] <= -fc.num] <- "Down"
  return(x)
}
