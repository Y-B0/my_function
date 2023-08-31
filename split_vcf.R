library(vcfR)
vcf<-read.vcfR("file/lnd.filtered.final.vcf")
split_vcf<-function(vcf,out_dir){
  meta<-vcf@meta
  fix<-vcf@fix
  gt<-vcf@gt
  result<-list()
  for (i in 2:ncol(gt)) {
    result[[i-1]]<-vcf
    result[[i-1]]@meta<-meta
    result[[i-1]]@fix<-fix
    result[[i-1]]@gt<-gt[,c(1,i)]
    write.vcf(result[[i-1]],file = paste(out_dir,colnames(gt)[i],".vcf.gz",sep = ""))
  }
}
split_vcf(vcf,"./file/vcf/")

