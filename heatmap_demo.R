##heatmap
heatmap_plot<-function(x,exprs,fc.name="logFC",nhead=50,ntail=50){
  #x:the datafram of deg product from upstream deg object
  #y:the normlized matrix of result from upstream deg object
  #sig:the
  library(dplyr)
  library(ComplexHeatmap)
  ##
  x<-x[order(x[,fc.name],decreasing = T),]
  g_name<-rownames(rbind(head(x,nhead),tail(x,ntail)))
  exprs<-exprs[rownames(exprs)%in%intersect(rownames(exprs),g_name),]
  exprs<-scale(exprs,center = T,scale = T)
  Heatmap(exprs,row_dend_reorder = TRUE,row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),cluster_columns=F)
}