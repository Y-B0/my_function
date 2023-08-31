library(stringr)
library(biomaRt)
library(Hmisc)
library(psych)
library(org.Hs.eg.db)
library(BiocGenerics)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(RColorBrewer)

deg <- read.delim("F:/deg_bye_0718/deg.txt")
deg<-deg[,c(1,4)]
deg$gene_name<-str_to_upper(deg$gene_name)


mart<-useDataset("hsapiens_gene_ensembl",useMart("ensembl"))
symbol<-getBM(attributes = c("entrezgene_id","hgnc_symbol"),
              filters = "hgnc_symbol",values = deg$gene_name,mart=mart)
symbol<-symbol[!duplicated(symbol$hgnc_symbol),]
tmp<-merge(deg,symbol,by.x="gene_name",by.y="hgnc_symbol",all=T)
tmp<-tmp[order(tmp$log2FoldChange,decreasing = T),]

geneList <- tmp[,2]
names(geneList) = as.character(tmp[,3])
geneList<-sort(geneList,decreasing = T)

Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
R.utils::setOption("clusterProfiler.download.method","auto")
KEGG_gseresult <- kegggsea<-gseKEGG(geneList,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH",
                                    verbose = T,seed = T)
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)

write.table (Go_gseresult, file ="Go_gseresult.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="KEGG_gseresult.csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="Go_Reactomeresult.csv", sep =",", row.names =TRUE)

cp<-Go_gseresult
Go_gseresult@result<-Go_gseresult@result[order(Go_gseresult@result$NES),]

Go_gseresult@result<-Go_gseresult@result[Go_gseresult@result$ID=="GO:0090257",]
gseaplot2(Go_gseresult,1,base_size = 15,color = pal_d3()(1),pvalue_table = T) 

gp(Go_gseresult,1,base_size = 15,color = pal_d3()(1),pvalue_table = T) 
