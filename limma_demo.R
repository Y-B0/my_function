##limma 
limma_after_normalized<-function(exp,group,compare,normalize=F,log2=F,merge=F,rna.count=F,symbol.name="SYMBOL"){
  ## exp is expression data (if merge==T,exp need a col to contain gene symbol)
  ## group can be multi group ,and first col is sample name; second col is condition; coef, and multi used to select multi group
  ## compare used to identify which group vs which group (usually <conventional> - <control>). also use for multi group condition.
  
  if (merge==T) {
    expr = aggregate(expr, by=list(expr[,symbol.name]),mean)
    expr = expr[!is.na(expr[,symbol.name]) & expr[,symbol.name]!="",]
    row.names(expr) = expr[,symbol.name]
    expr<-expr[,-grep(symbol.name,colnames(expr))]
  }
  
  if (rna.count==T) {
    expr<- voom(exprSet,design,normalize="quantile")$E
  }
  
  if (normalize==T) {
    expr<-normalizeBetweenArrays(expr)
  }
  if (log2==T) {
    expr<-log2(expr+1)
  }
  
  design <- model.matrix(~0+factor(group[,2]))
  colnames(design) <- unique(group[,2])
  rownames(design) <- group[,1]
  exp<-exp[,rownames(design)]

  fit <- lmFit(exp,design)
  contrast.matrix <- makeContrasts(compare,levels = design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  output <- topTable(fit2,coef = compare,n=Inf)
  return(output)
}
