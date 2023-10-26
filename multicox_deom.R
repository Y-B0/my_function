multicox<-function(expr,covariates){
  library(survival)
  ### x is the data frame contain the expression of variates and time, status. keep in mind ,the colname of time and status must be 'time' and 'status'
  ### covariates are the variates that needed to analysis
  res.cox <- coxph(eval(parse(text = paste("Surv(time, status) ~ ",paste(covariates,collapse = " + "),""))), data =  expr)
  x <- summary(res.cox)
  
  pvalue=signif(as.matrix(x$coefficients)[,5],2)
  HR=signif(as.matrix(x$coefficients)[,2],2)
  low=signif(x$conf.int[,3],2)
  high=signif(x$conf.int[,4],2)
  
  multi_res=data.frame(p.value=pvalue,
                       HR=paste(HR," (",low,"-",high,")",sep=""),
                       stringsAsFactors = F)
  return(multi_res)
  
}
