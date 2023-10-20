my_custom_theme <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.background=element_rect(fill = NA)
    )
}