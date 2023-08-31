# BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## 下面是定义基因长度为 非冗余exon长度之和
if(!file.exists('gene_length.Rdata')){
  exon_txdb=exons(txdb)
  genes_txdb=genes(txdb)
  
  o = findOverlaps(exon_txdb,genes_txdb)
  o
  t1=exon_txdb[queryHits(o)]
  t2=genes_txdb[subjectHits(o)]
  t1=as.data.frame(t1)
  t1$geneid=mcols(t2)[,1]
  
  # 如果觉得速度不够，就参考R语言实现并行计算
  # http://www.bio-info-trainee.com/956.html
  g_l = lapply(split(t1,t1$geneid),function(x){
    # x=split(t1,t1$geneid)[[1]]
    head(x)
    tmp=apply(x,1,function(y){
      y[2]:y[3]
    })
    length(unique(unlist(tmp)))
  })
  head(g_l)
  g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))
  
  save(g_l,file = 'gene_length.Rdata')
} 

load(file = 'gene_length.Rdata')
head(g_l)