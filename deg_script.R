library(readxl)
library(limma)
library(magrittr)
library(stringr)
library(edgeR)


############################
#exp need colname for sample rownames for gene
#group need control before case
#if count alread produce 2 log2, need return 2 count first, use return.count=T
#
############################

##limma used for alread normalized 2 deg
limma_after_normalized<-function(exp,group,return.count=F){
  if (return.count==F) {
    exp<-2^exp-1
  }
  design <- model.matrix(~0+factor(group))
  colnames(design) <- unique(group)
  rownames(design) <- colnames(exp)
  fit <- lmFit(exp,design)
  fit2 <- eBayes(fit)
  output <- topTable(fit2, coef=2,n=Inf)
  return(output)
}


##edgeR from count to deg
edgeR_count2deg<-function(exp,group,return.count=F){
  if (return.count==T) {
    exp<-2^exp-1
  }
  y<-DGEList(counts = exp,group = group)
  keep<-filterByExpr(y)
  y<-y[keep,,keep.lib.size=F]
  y<-calcNormFactors(y)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  y<-estimateDisp(y)
  et<-exactTest(y)
  et<-topTags(et,n=Inf)
  et<-as.data.frame(et)
  et<-cbind(rownames(et),et)
  colnames(et)<-c("id","log2FC","log2CPM","Pvalue","adj.P.Val")
  return(list(et,logcpm))
}

##limma used for count 2 deg
limma_count2deg<-function(exp,group,return.count=F){
  if (return.count==T) {
    exp<-2^exp-1
  }
  dge <- DGEList(counts = exp)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  design <- model.matrix(~group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(exp)
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend=TRUE)
  output <- topTable(fit, coef=2,n=Inf)
  colnames(output)<-c("log2FC","AveExpr","t","Pvalue","adj.P.Val","B")
  output$id<-rownames(output)
  return(list(output,logCPM))
}

#downstream analysis was based on edgeR_count2deg, other function need to keep colname as same.
##annotation


ensembl2symbol<-function(deg,id="id",from="ensembl_gene_id",to="hgnc_symbol"){
  library(biomaRt)
  ensembl=useMart("ensembl")
  hsa_ensembl = useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)
  gene_symbol<-getBM(attributes = c(from,to),filters = c(from),values = deg[,id],mart = hsa_ensembl)
  deg<-merge(deg,gene_symbol,by.x=id,by.y=from)
  return(deg)
}

##volcano plot
#### base on the regult of edgeR, default p<0.01 and abs(log2fc)>2
volcano<-function(x,adj="P.Value",log="logFC"){
  library(ggplot2)
  x$sig[(x[,adj] > 0.05|x[,adj]=="NA")|(x[,log] < 0.585)& x[,log] > -0.585] <- "No"
  x$sig[x[,adj] <= 0.05 & x[,log] >= 0.585] <- "Up"
  x$sig[x[,adj] <= 0.05 & x[,log] <= -0.585] <- "Down"
  ggplot()+theme_bw()+xlim(-10,10)+
    geom_point(aes(x=x[,log],y=-1*log10(x[,adj]),color=x$sig))+theme(text = element_text(size=20))+
    labs(x="log2(FoldChange)",y="-log10(FDR)")+scale_color_manual(name="",values =c("#0072B5","grey","#BC3C28"))
}

##heatmap
heatmap_plot<-function(x,logcpm,adj="Pvalue",nhead=100,ntail=100){
  #x:the datafram of deg product from upstream deg object
  #y:the normlized matrix of result from upstream deg object
  #make sure the rownames logcpm equl to x$id
  library(dplyr)
  library(ComplexHeatmap)
  et<-x[abs(x$log2FC)>=1&x$`adj`<0.01,] %>%.[order(.$log2FC),]
  et<-rbind(head(et,nhead),tail(et,ntail))%>%.[order(.$log2FC,decreasing = T),]
  et<-et[et$hgnc_symbol!="",]
  logcpm<-logcpm[rownames(logcpm)%in%intersect(rownames(logcpm),et$id),]
  rownames(logcpm)<-et$hgnc_symbol
  Heatmap(logcpm,row_dend_reorder = TRUE,row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),cluster_columns=T)
}

##enrichment
enrich<-function(gene_symbol,pvalue=0.05,n=10){
  ##GO??????
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(Hmisc)
  library(ggsci)
  gene.df <- bitr(gene_symbol, fromType="SYMBOL",
                  toType="ENTREZID",
                  OrgDb = "org.Hs.eg.db")
  go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
  write.table(go,file = "GO_ALL.txt",sep = "\t",quote = F,row.names = F)
  go<-go@result
  go<-split(go[go$pvalue<pvalue,],go$ONTOLOGY)%>%lapply(.,head,n=n)%>%Reduce(rbind,.)
  go<-go[order(go$Count,decreasing = T),]
  go<-go[order(go$ONTOLOGY,decreasing = F),]
  go$Description <- factor(go$Description,levels = rev(go$Description))
  ggplot(data = go, 
         aes(x = Description, y = Count,fill = ONTOLOGY))+ 
    geom_bar(stat = "identity",width = 0.9)+ 
    coord_flip()+theme_bw()+ 
    scale_x_discrete(labels = function(x) str_wrap(x,width = 80))+ 
    labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ 
    theme(axis.title = element_text(size = 13), 
          axis.text = element_text(size = 11),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))+
    scale_fill_aaas()
  ggsave("go.pdf",plot = p,width = 8,height = 10)
  ##KEGG

  kegg <- enrichKEGG(gene.df$ENTREZID, keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE)
  kegg=DOSE::setReadable(kegg, OrgDb='org.Hs.eg.db', keyType = 'ENTREZID')
  kegg$GeneRatio<-lapply(kegg$GeneRatio, function(x){eval(parse(text=x))})%>%unlist()
  write.table(data.frame(kegg),"KEGG_ALL.txt",sep = "\t",quote = F,row.names = F)
  ggplot(kegg[1:n,], aes(x=GeneRatio, y=Description,size=Count, color=pvalue)) + geom_point() + 
    scale_colour_gradient(low="green",high="red") + labs(color=expression(padj),size="Gene number", x="GeneRatio",y="Pathway name",title="KEGG Pathway enrichment")
  ggsave("kegg.pdf",plot = p,width = 8,height = 6)
  ##cytoband
  cytoband_gmt <- read.gmt("D:/project/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c1.all.v7.4.entrez.gmt")
  cytoband<-enricher(gene.df$ENTREZID,TERM2GENE = cytoband_gmt)
  write.table(as.data.frame(cytoband@result),file = "cytoband.txt",row.names = F,col.names = T,quote = F,sep = "\t")
  
  return(list(GO=go@result,KEGG=kegg@result,Cytoband=cytoband@result))
}

##gsea

gsea<-function(deg){
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(Hmisc)
  gene=bitr(deg$hgnc_symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  data_all <- merge(deg,gene,by.x="hgnc_symbol",by.y="SYMBOL")%>% arrange(desc(log2FC))
  geneList = data_all$log2FC
  names(geneList) <- data_all$ENTREZID
  
  gse.GO <- gseGO(
    geneList, 
    ont = "BP",  
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  write.csv(gse.GO[[1]],"gse_GO.csv",row.names = F)
  gse.KEGG <- gseKEGG(geneList, 
                      organism = "hsa", 
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH")
  return(list(gse.GO=gse.GO,gse.KEGG=gse.KEGG))
  write.csv(gse.KEGG[[1]],"gse_KEGG.csv",row.names = F)
}



