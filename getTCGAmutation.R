

getTCGAmutation <-function(project){
  library(TCGAmutations)
  
  try(dir.create("output_mute"))
  
  mutation<-TCGAmutations::tcga_load(study = strsplit(project,split = "-")[[1]][2], source = "Firehose")
  save(mutation,paste0("output_mute/",project,"_mute.rdata")) 
  
}