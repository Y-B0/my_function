#########################################################
#  limma used for arrays
#  Deseq2 and limma voom used for RNA-seq need raw count 
# 
#########################################################

library(GEOquery)
library(limma)
GSE="GSE62402"
gset <- getGEO(GSE, destdir=".",
               AnnotGPL = F,     ## annotation file
               getGPL = F)       ## plamfor file
save.image(paste(GSE,".RData",sep = ""))
## read count
exprSet <- exprs(gset[[1]])
########################################## 
#boxplot(exprSet)
#exprSet<-normalizeBetweenArrays(exprSet) #option. make different group has the same distribution
#exprSet<-log2(exprSet+1) #options
##########################################
## group data
pData <- pData(gset[[1]])
## probe data
gpl <- getGEO(gset[[1]]@annotation, destdir = "./") ## or gpl<-getIDs(gset[[1]]@annotation)
gpl_anno= Table(gpl)
ID2Symbol = data.frame("ID"=gpl_anno[,1],"SYMBOL"=gpl_anno[,2])
## expression process
exprSet_temp = as.data.frame(exprSet)
exprSet_temp$ID = row.names(exprSet_temp)
exprSet_temp = merge(exprSet_temp, ID2Symbol, by.x="ID", by.y="ID", all.x=T)
exprSet_temp = exprSet_temp[,-1]
exprSet_temp = aggregate(exprSet_temp, by=list(exprSet_temp[,"SYMBOL"]),max)
exprSet_temp = exprSet_temp[!is.na(exprSet_temp$SYMBOL) & exprSet_temp$SYMBOL!="",]
row.names(exprSet_temp) = exprSet_temp$SYMBOL
exprSet_temp = exprSet_temp[,2:(ncol(exprSet_temp)-1)]
exprSet = exprSet_temp

## group 
sample <- pData$geo_accession[grep("lacZ.*|Mef2c.*",pData$title)]
group<-pData$title[grep("lacZ.*|Mef2c.*",pData$title)]
group<-gsub(" | knockdown.*","",group)
sample_group <- data.frame(sample, group)
design <- model.matrix(~0+factor(sample_group$group))
colnames(design) <- levels(factor(sample_group$group))
rownames(design) <- sample_group$sample
exprSet<-exprSet[,rownames(design)]

contrast.matrix = makeContrasts(eval(paste(colnames(design),collapse = "-")),levels = design)


## deg
fit = lmFit(exprSet,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
x = topTable(fit2,coef = 1,n=Inf, adjust.method = "BH", sort.by="P")

## annotation
ensembl<-useMart("ensembl")
hsa_ensembl = useDataset(dataset="mmusculus_gene_ensembl", mart=ensembl)
gene_symbol<-getBM(attributes = c("entrezgene_id","ensembl_gene_id"),filters = "entrezgene_id",values = rownames(x),mart = hsa_ensembl)
x<-x[match(gene_symbol$entrezgene_id,rownames(x)),]
x$ensg<-gene_symbol$ensembl_gene_id

col<-"adj.P.Val"
fc<-"logFC"
x$sig[(x[,col] > 0.05|x[,col]=="NA")|(x[,fc] < 0.585)& x[,fc] > -0.585] <- "no"
x$sig[x[,col] <= 0.05 & x[,fc] >= 0.585] <- "up"
x$sig[x[,col] <= 0.05 & x[,fc] <= -0.585] <- "down"



