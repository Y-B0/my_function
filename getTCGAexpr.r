


getTCGAexpr <- function(project){
  library(TCGAbiolinks)
  
  try(dir.create("output_expr"))

  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, files.per.chunk = 100)

  GDCprepare(query,save = T,save.filename = paste0("output_expr/",project,"_expr.rdata"))
  
  load(file = paste0("output_expr/",project,"_expr.rdata"))
  
  se <- data
  
  clin_info <- as.data.frame(colData(se))
  save(clin_info, file = paste0("output_expr/",project, "_clinical.rdata"))
  
  rowdata <- rowData(se)
  
  se_mrna <- se[rowdata$gene_type == "protein_coding",]
  se_lnc <- se[rowdata$gene_type == "lncRNA",]
  
  
  expr_counts_mrna <- assay(se_mrna,"unstranded")

  expr_tpm_mrna <- assay(se_mrna,"tpm_unstrand")

  expr_fpkm_mrna <- assay(se_mrna,"fpkm_unstrand")
  
  expr_counts_lnc <- assay(se_lnc,"unstranded")
  
  expr_tpm_lnc <- assay(se_lnc,"tpm_unstrand")
  
  expr_fpkm_lnc <- assay(se_lnc,"fpkm_unstrand")
  
  symbol_mrna <- rowData(se_mrna)$gene_name
  
  symbol_lnc <- rowData(se_lnc)$gene_name
  

  mrna_expr_counts <- cbind(data.frame(symbol_mrna),as.data.frame(expr_counts_mrna))
  mrna_expr_counts <- mrna_expr_counts %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  save(mrna_expr_counts, file = paste0("output_expr/",project, "_mrna_expr_counts.rdata"))
  
  
  mrna_expr_tpm <- cbind(data.frame(symbol_mrna),as.data.frame(expr_tpm_mrna))
  mrna_expr_tpm <- mrna_expr_tpm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  save(mrna_expr_tpm, file = paste0("output_expr/",project, "_mrna_expr_tpm.rdata"))
  
  
  mrna_expr_fpkm <- cbind(data.frame(symbol_mrna),as.data.frame(expr_fpkm_mrna))
  mrna_expr_fpkm <- mrna_expr_fpkm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  save(mrna_expr_fpkm, file = paste0("output_expr/",project, "_mrna_expr_fpkm.rdata"))
  
  
  lncrna_expr_counts <- cbind(data.frame(symbol_lnc),as.data.frame(expr_counts_lnc))
  lncrna_expr_counts <- lncrna_expr_counts %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_lnc,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_lnc") %>% 
    as.data.frame()
  save(lncrna_expr_counts, file = paste0("output_expr/",project, "_lncrna_expr_counts.rdata"))
  
  
  lncrna_expr_tpm <- cbind(data.frame(symbol_lnc),as.data.frame(expr_tpm_lnc))
  lncrna_expr_tpm <- lncrna_expr_tpm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_lnc,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_lnc") %>% 
    as.data.frame()
  save(lncrna_expr_tpm, file = paste0("output_expr/",project, "_lncrna_expr_tpm.rdata"))
  
  
  lncrna_expr_fpkm <- cbind(data.frame(symbol_lnc),as.data.frame(expr_fpkm_lnc))
  lncrna_expr_fpkm <- lncrna_expr_fpkm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_lnc,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_lnc") %>% 
    as.data.frame()
  save(lncrna_expr_fpkm, file = paste0("output_expr/",project, "_lncrna_expr_fpkm.rdata"))
  
  
}













