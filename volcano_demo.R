volcano_plot<-function(x, p.name = "P.Value", fc.name = "logFC", p.value = 0.05, fc.value = 0.585,file.name=NULL,plot.name=NULL) {
  library(ggplot2)

  x$sig[(x[, p.name] > p.value | x[, p.name] == "NA") | (x[, fc.name] < fc.value) & x[, fc.name] > -fc.value] <- "No"
  x$sig[x[, p.name] <= p.value & x[, fc.name] >= fc.value] <- "Up"
  x$sig[x[, p.name] <= p.value & x[, fc.name] <= -fc.value] <- "Down"

  p <- ggplot() +
    theme_bw() +
    xlim(-10, 10) +
    geom_point(aes(x = x[, fc.name], y = -1 * log10(x[, p.name]), color = x$sig)) +
    theme(text = element_text(size = 20)) +
    labs(x = "log2(FoldChange)", y = paste("-log10(", p.name, ")", sep = "")) +
    scale_color_manual(name = "", values = c("#0072B5", "grey", "#BC3C28"))


  if (!is.null(file.name)) {
    write.table(data.frame(Symbol = rownames(x), x), file = file.name, sep = "\t", quote = F, row.names = F, col.names = T)
  }

  if (!is.null(plot.name)) {
    ggsave(plot.name, p, "pdf", width = 7, height = 4.5)
  }

  return(list(data = x, plot = p))
}
