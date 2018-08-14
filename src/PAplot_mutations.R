# Presence/Absence plot for mutations
PAplot_mutations <- function(df) {
  cols <- c("1" = "#95cbee","0" = "#c9e2f6")
  
  p <- ggplot(df, aes(gene, ref, fill = result)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = cols,
                      labels = c("Mutated","Wild Type"),
                      breaks = c("1","0"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 18),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank())+
    coord_fixed(0.5)
  
  ggsave("mutation_plot.tiff",
         p,
         device = "tiff",
         dpi = 300,
         width = 40,
         height = 60,
         units = "cm")
}