wgcna_deom<-function(exprs,feature,type="signed",corType = "bicor"){
  ## exprs: matrix of expression, row means gene and col means samples
  ## feature: dataframe of tarit data, row means sample and col means features
  library(magrittr)
  library(WGCNA)
  library(stringr)
  
  goo.gene <- goodSamplesGenes(t(exprs), verbose = 3);
  if (!goo.gene$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!goo.gene$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(exprs)[!goo.gene$goodGenes], collapse = ", ")));
    if (sum(!goo.gene$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(exprs)[!goo.gene$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    exprs = t(exprs)[goo.gene$goodSamples, goo.gene$goodGenes]
  }
  
  exprs<-t(exprs)
  
  m.mad <- apply(exprs,1,mad)
  dataExpr <- exprs[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]%>%t(.)
  
  powers<-1:20
  nSamples=nrow(dataExpr)
  ngene <- ncol(dataExpr)
  type="signed"
  corType = "bicor"
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  robustY = ifelse(corType=="pearson",T,F)
  sft <- pickSoftThreshold(dataExpr, powerVector = powers)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], 
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", 
       ylab = "SFT, signed R^2", type = "n", 
       main = paste("Scale independence"))
  text(sft$fitIndices[, 1], 
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = 1:20, col = "red")
  abline(h = 0.85, col = "red")
  plot(sft$fitIndices[, 1], 
       sft$fitIndices[, 5], 
       type = "n", 
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity", 
       main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], 
       sft$fitIndices[, 5], 
       labels = powers, col = "red")
  power_plot<-recordPlot()
  power<-sft[[1]]
  
  sampleTree = hclust(dist(dataExpr), method = "average")
  
  net = blockwiseModules(dataExpr, power = power, maxBlockSize = 5000,
                         TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.35,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, corType = corType, 
                         loadTOMs=TRUE,maxPOutliers=maxPOutliers,
                         saveTOMFileBase = paste0("result", ".tom"),
                         verbose = 3)
  
  ##plot cluster
  
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  cluster_plot<-plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,)
  
  
  MEs = net$MEs
  MEs_col = net$MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  
  module_heatmap<-plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                        xLabelsAngle = 90,)
  
  sampleName = rownames(dataExpr)
  traitData = feature[match(sampleName, rownames(feature)), ]
  
  if (corType=="pearsoon") {
    modTraitCor = cor(MEs_col, traitData, use = "p")
    modTraitP = corPvalueStudent(modTraitCor, nSamples)
  } else {
    modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
    modTraitCor = modTraitCorP$bicor
    modTraitP   = modTraitCorP$p
  }
  textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
  dim(textMatrix) = dim(modTraitCor)
  
  tarit_heatmap<-labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
                 yLabels = colnames(MEs_col), 
                 cex.lab = 0.5, 
                 xLabelsAngle = 0,
                 ySymbols = colnames(MEs_col), colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 setStdMargins = FALSE, textMatrix = textMatrix,
                 cex.text = 0.5, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  if (corType=="pearsoon") {
    geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(
      as.matrix(geneModuleMembership), nSamples))
  } else {
    geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
    geneModuleMembership = geneModuleMembershipA$bicor
    MMPvalue   = geneModuleMembershipA$p
  }
  
  
  if (corType=="pearsoon") {
    geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
    geneTraitP = as.data.frame(corPvalueStudent(
      as.matrix(geneTraitCor), nSamples))
  } else {
    geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
    geneTraitCor = as.data.frame(geneTraitCorA$bicor)
    geneTraitP   = as.data.frame(geneTraitCorA$p)
  }
  
  ##plot MM-GS
  modNames = substring(names(MEs_col), 3)
  dir.create("./GM_MM/")
  setwd("./GM_MM/")
  for (i in modNames) {
    module = i
    for (j in colnames(traitData)) {
      pheno=j
      module_column = match(module, modNames)
      pheno_column = match(pheno,colnames(traitData))
      moduleGenes = moduleColors==module
      pdf(file = paste(module,pheno,".pdf",sep = "_"),width = 5,height = 5)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                         abs(geneTraitCor[moduleGenes, pheno_column]),
                         xlab = paste("Module Membership in", module, "module and",pheno),
                         ylab = "Gene significance for body weight",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
  setwd("../")
  
  gene<-as.data.frame(net$colors)
  gene$`net$colors`<-paste("ME",gene$`net$colors`)
  gene$color<-paste0("ME",labels2colors(as.numeric(str_replace_all(gene$`net$colors`,"ME",""))))
  gene$gene<-names(net$colors)
  
  module_g<-gene[,2:3]
  module_g$color<-factor(module_g$color,level = names(table(module_g$color)[order(table(module_g$color))]))
  module_g<-module_g[order(module_g$color),]
  module_g<-split(module_g,module_g$color)
  module_g<-lapply(module_g,function(x){x[,2]})
  module_g<-do.call(cbind, lapply(lapply(module_g, unlist), `length<-`, max(lengths(module_g))))
  module_g[is.na(module_g)]<-""
  write.table(module_g,file = "module_gene.txt",sep = "\t",quote = F,row.names = F,col.names = T)
  
  
  load(net$TOMFiles[1], verbose=T)
  TOM = TOMsimilarityFromExpr(dataExpr, power = power);
  TOM <- as.matrix(TOM)
  dissTOM = 1-TOM
  plotTOM = dissTOM^7
  diag(plotTOM) = NA
  #TOMplot(plotTOM, net$dendrograms, moduleColors, 
  #        main = "Network heatmap plot, all genes")
  dir.create("./network")
  for (i in modNames) {
    module=i
    probes = colnames(dataExpr)
    inModule = is.finite(match(moduleColors, module))
    modProbes = probes[inModule]
    modTOM = TOM[inModule, inModule]
    dimnames(modTOM) = list(modProbes, modProbes)
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = paste("./network/",module,"model", ".edges.txt", sep=""),
                                   nodeFile = paste("./network/",module,"model", ".nodes.txt", sep=""),
                                   weighted = TRUE, threshold = 0,
                                   nodeNames = modProbes, nodeAttr = moduleColors[inModule])
  }
  save.image("wgcna.RData")
  
}