unicox<-function(expr,covariates){
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
                           HR.merge <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(HR,HR.confint.lower,HR.confint.upper,HR.merge,p.value)
                           names(res)<-c("HR","HR.lower","HR.upper","HR.merge","p.value")
                           return(res)
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res<-data.frame(Module=covariates,res)
  return(res)
}

forestplot(res[,c(1,2,5,6)],
           mean=res$HR,   #告诉函数，表格第2列为HR，它要变成森林图的小方块
           lower=res$HR.confint.lower,  #告诉函数表格第3列为5%CI，
           upper=res$HR.confint.upper,  #表格第5列为95%CI，它俩要化作线段，穿过方块
           zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
           boxsize=0.2,       #设置小黑块的大小
           graph.pos=2,       #森林图应插在图形第2列
           #graph.pos="right",       #森林图应插在图形最右边
           graphwidth = unit(.3,"npc"), # 图的宽度
           xticks=c(0.5,1,2,3,4,5,6,7,8),
           fontsize = 16,
           txt_gp=fpTxtGp(
             label=gpar(cex=1),
             ticks=gpar(cex=1), 
             xlab=gpar(cex=1), 
             title=gpar(cex=1.5))
)
