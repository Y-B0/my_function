merge_TCGA <- function(manifest, json, path,gene_type="protein_coding"){
  library(rjson)
  json <- jsonlite::fromJSON(txt = json)
  id <- rep(0, length(json))
  sample_id <- rep(0, length(json))
  manifest <- read.delim(manifest)
  for (i in 1:nrow(json)) {
    id[i] = json$file_id[i]
    sample_id[i] = json$associated_entities[[i]]$entity_submitter_id
  }
  sample_matrix = data.frame(id = id, sample_id = sample_id)
  manifest<-manifest[match(sample_matrix$id,manifest$id),]
  x_merge<-NULL
  exp<-apply(manifest,1,function(x){
    filename <- paste(path,x[1],x[2],sep = "/")
    x <- read.delim(filename, comment.char="#")%>%.[-(1:4),]
    x <- x[!duplicated(x$gene_name),]
    rownames(x)<-x[,"gene_name"]
    x<-x[x$gene_type==gene_type,c("gene_name","tpm_unstranded")]
    x
  })
  exp<-Reduce(cbind,exp)
  exp<-exp[,colnames(exp)!="gene_name"]
  write.table(exp,file = "exp_count.txt",quote = F,sep = "\t",row.names = T,col.names = T)
}
