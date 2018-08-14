# Presence/absence plot for acquired genes
PAplot_acquired_genes <- function(df) {
  cols <- c("1" = "#95cbee","0" = "#c9e2f6")
  
  p <- ggplot(df, aes(gene, ref, fill = result)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = cols,
                      labels = c("Present","Absent"),
                      breaks = c("1","0"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank())+
    coord_fixed()
  
  ggsave("acquired_genes_plot.tiff",
         p,
         device = "tiff",
         dpi = 100,
         width = 20,
         height = 20,
         units = "cm")
}