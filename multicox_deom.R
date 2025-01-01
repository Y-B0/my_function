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
  HR.merge <- paste0(HR, " (", low, "-", high, ")")
  multi_res=data.frame(
    Model=covariates,
    HR=HR,
    HR.lower=low,
    HR.high=high,
    HR.merge=HR.merge,
    p.value=pvalue,
    stringsAsFactors = F)
  return(multi_res)

  
}

forestplot(multi_res[,c(1,2,5,6)],
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
