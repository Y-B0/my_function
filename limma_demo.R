limma_after_normalized<-function(exp,group,compared,normalize=F,log2=F,merge=F,rna.count=F,symbol.name="SYMBOL"){
  ## exp is expression data (if merge==T,exp need a col to contain gene symbol)
  ## group can be multi group ,and first col is sample name; second col is condition; coef, and multi used to select multi group
  ## compared used to identify which group vs which group (usually <conventional> - <control>). also use for multi group condition.
  print(compared)
  if (merge==T) {
    exp = aggregate(exp, by=list(expr[,symbol.name]),mean)
    exp = expr[!is.na(exp[,symbol.name]) & exp[,symbol.name]!="",]
    row.names(expr) = exp[,symbol.name]
    exp<-exp[,-grep(symbol.name,colnames(exp))]
  }
  
  if (rna.count==T) {
    design <- model.matrix(~0 + factor(group[, 2]))
    exp<- voom(exp,design,normalize="quantile")$E
  }
  
  if (normalize==T) {
    exp<-normalizeBetweenArrays(exp)
  }
  if (log2==T) {
    exp<-log2(exp+1)
  }
  
  design <- model.matrix(~0 + factor(group[, 2]))
  colnames(design) <- unique(group[, 2])
  rownames(design) <- group[, 1]
  exp <- exp[, rownames(design)]
  
  fit <- lmFit(exp, design)
  contrast.matrix <- makeContrasts(contrasts=compared, levels = colnames(design))  # 修改这里
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  output <- topTable(fit2, coef = compared, n = Inf)
  return(cbind(exp,output))
}
