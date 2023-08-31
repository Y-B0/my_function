#pcatest
load("basic_data.RData")
#test<-PCA(BRCA.mrna.select[,c(na.omit(gene[,2]))])
#tmp<-get_pca_ind(test)
#tmp1<-as.data.frame(tmp[[2]])
for (i in 1:ncol(gene)) {
  data<-BRCA.mrna.select[,c(na.omit(gene[,i]))]
  test<-prcomp(data,center = T)
  tmp<-as.data.frame(test$x)
  ##alter symbol based on samle
  #ifelse(sum(tmp[,1]>0)>520,tmp[,1]<-abs(tmp[,1]),tmp[,1]<-(-abs(tmp[,1])))
  ##alter symbol based on gene
  tmp[,1]<-ifelse(rowSums(data>0)>(ncol(data)/2),abs(tmp[,1]),(-abs(tmp[,1])))
  pca[,i]<-tmp[,1]
  
}
#save(pca,file = "pca_prcomp_no_center.RData")





##pca heatmap plot
i<-8
test<-prcomp(BRCA.mrna.select[,c(na.omit(gene[,i]))],center = T)
tmp<-as.data.frame(test$x)
#tmp1<-svd(BRCA.mrna.select[,c(na.omit(gene[,2]))])%>%.[[2]]
data<-cbind(BRCA.mrna.select[,c(na.omit(gene[,i]))],tmp[,1],pca[,i])
data<-apply(data, 2,as.numeric)
colnames(data)[(ncol(data)-1):ncol(data)]<-c("pca","svd")
Heatmap(as.matrix(data),show_row_names = F,show_column_names = F,cluster_columns = F,cluster_rows = F)


    


  
