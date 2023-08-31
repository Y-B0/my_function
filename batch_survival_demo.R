## survival
library(survival)
library(survminer)
####### data read and sample fileter
load("READ.RData")
rm(list = ls()[!ls()%in%c("logcpm","group","deg")])
rownames(logcpm)<-str_split(rownames(logcpm),"\\.",simplify = T)[,1]
deg<-deg[deg$hgnc_symbol!="",]
deg<-deg[!duplicated(deg$hgnc_symbol),]
deg<-deg[!duplicated(deg$id),]
deg<-deg[deg$id%in%intersect(rownames(logcpm),deg$id),]
logcpm<-logcpm[rownames(logcpm)%in%intersect(rownames(logcpm),deg$id),]
rownames(logcpm)<-deg$hgnc_symbol
TFs <- read.delim("Integrated_meanRank.tsv")%>%head(.,n=10)
logcpm<-logcpm[intersect(TFs$TF,rownames(logcpm)),group$sample[group$group=="T"]]
logcpm<-t(logcpm)
TCGA.READ.survival <- read.delim("~/Desktop/project/TCGA/TCGA-READ.survival.tsv")
TCGA.READ.survival$sample<-gsub("-",".",TCGA.READ.survival$sample)
logcpm<-logcpm[intersect(rownames(logcpm),TCGA.READ.survival$sample),]
rownames(TCGA.READ.survival)<-TCGA.READ.survival$sample
TCGA.READ.survival<-TCGA.READ.survival[intersect(rownames(logcpm),TCGA.READ.survival$sample),]
#######
####### merge data and survival data
data<-cbind(TCGA.READ.survival[,c(2,4)],logcpm)
colnames(data)[1:2]<-c("status","time")
res.cut <- surv_cutpoint(data,
                         time = "time",
                         event = "status",
                         variables = colnames(data)[3:ncol(data)] #select gene
)
surv_group <- surv_categorize(res.cut)
surv_group<-melt(surv_group,id.vars=c("status","time"))
fit<- survfit(Surv(time, status) ~ value, data = surv_group)
splots<-list()
for (i in colnames(surv_group)[3:ncol(surv_group)]){
  formu<-paste("Surv(time,status)~",i,sep = "")
  km_fit <- survfit(eval(parse(text=formu)), data=surv_group)
  splots[[i]]<-ggsurvplot(km_fit,
                          xlab = "Time_days",
                          ylab=i,##??????????????????
                          pval = T,
                          conf.int= F,
                          risk.table = T,
                          legend.title = i,
                          legend.labs = levels(surv_group[[i]]),##
                          surv.median.line = "hv",# ????????????
                          palette="lancet")
}
#######