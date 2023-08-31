
####DEP
{
  library("DEP")
  
  # Loading a package required for data handling
  library("dplyr")
  X1_Ori <- read_excel("~/Desktop/???????????????????????????????????????/HQ-BQ-ZQ20200619-FC3-TMT-XSZZ ???????????? ????????????/????????????/results/1_Statistics Analysis/1_ori.xlsx")
  X1_Ori <- X1_Ori[,c(1,6:14)]
  colnames(X1_Ori)[1]<-"ID"
  colnames(X1_Ori)[2:10]<-paste("X",colnames(X1_Ori)[2:10],sep = "")
  # Gene annotation
  library(org.Mm.eg.db)
  library(clusterProfiler)
  gene <- bitr(X1_Ori$ID, fromType = "UNIPROT", 
               toType = c("UNIPROT", "SYMBOL"), 
               OrgDb = org.Mm.eg.db)
  gene<-gene[match(intersect(gene$UNIPROT,X1_Ori$ID),gene$UNIPROT),]
  X1_Ori<-X1_Ori[match(intersect(gene$UNIPROT,X1_Ori$ID),X1_Ori$ID),]
  X1_Ori$name<-gene$SYMBOL
  # Are there any duplicated gene names?
  X1_Ori$name %>% duplicated() %>% any()
  
  # Make a table of duplicated gene names
  X1_Ori %>% group_by(name) %>% summarize(frequency = n()) %>% 
    arrange(desc(frequency)) %>% filter(frequency > 1)
  
  # Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
  X1_unique <- make_unique(X1_Ori, "name", "ID", delim = ";")
  
  # Are there any duplicated names?
  X1_unique$name %>% duplicated() %>% any()
  
  # Generate a SummarizedExperiment object using an experimental design
  group<-data.frame(label=colnames(X1_unique)[2:10],condition=rep(c("X4h","X6h","Xcontrol"),each=3),replicate=rep(3:1,times=3))
  X1_se <- make_se(X1_unique, 2:10, group)
  
  # Plot a barplot of the protein identification overlap between samples
  #plot_frequency(X1_se)
  
  # Filter for proteins that are identified in all replicates of at least one condition
  X1_filt <- filter_missval(X1_se, thr = 0)
  
  # Normalize the data
  X1_norm <- normalize_vsn(X1_filt) ##vsn normlize
  
  # Visualize normalization by boxplots for all samples before and after normalization
  plot_normalization(X1_filt, X1_norm)
  
  # Plot a heatmap of proteins with missing values
  plot_missval(X1_filt)
  
  # Plot intensity distributions and cumulative fraction of proteins with and without missing values
  plot_detect(X1_filt)
  
  
  # Test all possible comparisons of samples
  X1_diff <- test_diff(X1_norm, type = "all")
  dep <- add_rejections(X1_diff, alpha = 0.05, lfc = log2(1.5))
  
  
  # Plot the first and second principal components
  plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)
  
  library(ggplot2)
  library(ggrepel)
  
  data<-dep@assays@data@listData[[1]]%>%t(.)
  pca <- prcomp(data,center = F,scale. = F)
  PC<-as.data.frame(pca$x)
  
  summ <- summary(pca)
  xlab1 <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
  ylab1 <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  
  ggplot(data = PC,aes(x = PC1,y = PC2,color = group$condition))+
    stat_ellipse(aes(fill = group$condition),
                 type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
    geom_point(size = 3.5)+
    labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
    guides(fill = "none")+
    theme_bw()+
    scale_fill_manual(values = c("purple","orange","pink"))+
    scale_colour_manual(values = c("purple","orange","pink"))
  
  # Plot a heatmap of all significant proteins with the data centered per protein
  plot_heatmap(dep, type = "centered", kmeans = TRUE, 
               k = 6, col_limit = 4, show_row_names = T,
               indicate = c("condition", "replicate"))
  
  plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
               k = 6, col_limit = 10, show_row_names = FALSE)
  
  
  # Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
  plot_volcano(dep, contrast = "X4h_vs_X6h", label_size = 2, add_names = TRUE)
  plot_volcano(dep, contrast = "X4h_vs_Xcontrol", label_size = 2, add_names = TRUE)
  plot_volcano(dep, contrast = "X6h_vs_Xcontrol", label_size = 2, add_names = TRUE)
  # Plot a frequency plot of significant proteins for the different conditions
  plot_cond(dep)
  
  # Generate a results table
  data_results <- get_results(dep)
  
  # Number of significant proteins
  data_results %>% filter(significant) %>% nrow()
  
  # Generate a long data.frame
  df_long <- get_df_long(dep)
  save.image("dep.RData")
}

####limma

## Loading a package required for data handling
library("dplyr")
library(limma)
library(edgeR)
X1_Ori <- read_excel("~/Desktop/???????????????????????????????????????/HQ-BQ-ZQ20200619-FC3-TMT-XSZZ ???????????? ????????????/????????????/results/1_Statistics Analysis/1_ori.xlsx")
X1_Ori <- X1_Ori[,c(1,6:14)]
X1_Ori<-as.data.frame(X1_Ori)
rownames(X1_Ori)<-X1_Ori[,1]
X1_Ori<-X1_Ori[,-1]
X1_Ori<-X1_Ori[rowSums(X1_Ori>=1)>=ncol(X1_Ori)/2,]
X1_Ori<-X1_Ori[!duplicated(rownames(X1_Ori)),]

#X1_Ori<-log2(X1_Ori+1) ### only expresion in data,colname is sample ,rowname is gene

#
## Gene annotation
library(org.Mm.eg.db)
library(clusterProfiler)
gene <- bitr(rownames(X1_Ori), fromType = "UNIPROT", 
             toType = c("UNIPROT", "SYMBOL"), 
             OrgDb = org.Mm.eg.db)
gene<-gene[match(intersect(gene$UNIPROT,rownames(X1_Ori)),gene$UNIPROT),]
X1_Ori<-X1_Ori[match(intersect(gene$UNIPROT,rownames(X1_Ori)),rownames(X1_Ori)),]
raw_count<-X1_Ori
raw_count<-log2(raw_count+1)
group<-rep(c("4h","6h","control"),each=3)## single vector order same as sample 


##pca
pca <- prcomp(t(raw_count),center = T,scale = T)
PC<-as.data.frame(pca$x)

summ <- summary(pca)
xlab1 <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
group<-data.frame(label=colnames(X1_Ori),condition=rep(c("h4","h6","control"),each=3),replicate=rep(3:1,times=3))

ggplot(data = PC,aes(x = PC1,y = PC2,color = group$condition,shape=as.factor(group$replicate)))+
  stat_ellipse(aes(fill = group$condition),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",shape="Replicate",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))


##multi group DEG
group<-group$condition
design<-model.matrix(~0+factor(group))
colnames(design)<-levels(factor(unique(group)))
rownames(design)<-colnames(X1_Ori)

##
X1_Ori<-DGEList(counts = X1_Ori)
X1_Ori<-calcNormFactors(X1_Ori)
logCPM<-cpm(X1_Ori,log = T,prior.count = 3)
v<-voom(X1_Ori,design,normalize.method = "quantile")

deglist<-list()
for (x in lapply(combn(unique(group),2,simplify = F),paste,collapse = "-")) {
  contrast.matrix<-makeContrasts(x,levels=design)
  fit <- lmFit(v,design)
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  Output = topTable(fit2, coef=x, n=Inf)
  Output$gene<-gene$SYMBOL[match(rownames(Output),gene$UNIPROT)]
  Output$compare<-x
  deglist[[x]]<-Output
}
#save.image("deg.RData")

####vocanlno

