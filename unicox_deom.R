unicox <- function(expr,covariates){
  ### x is the data frame contain the expression of variates and time, status. keep in mind ,the colname of time and status must be 'time' and 'status'
  ### covariates are the variates that needed to analysis
  library(survival)
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time, status)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = expr)})
  
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           HR <-signif(x$coef[2], digits=2)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(p.value,HR)
                           names(res)<-c("p.value","HR")
                           return(res)
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  return(res)
}
