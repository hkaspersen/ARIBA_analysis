# Heatmap on number of mutations in each gene
create_mutation_heatmap <- function(df) {
  palette <- c("#e7f0fa","#c9e2f6","#95cbee",
               "#0099dc","#4ab04a","#ffd73e",
               "#eec73a","#e29421","#f05336",
               "#ce472e")
  
  p<- ggplot(df, aes(gene, ref, fill = sum)) +
    geom_tile(color = "white") +
    labs(fill = "Mutations") +
    scale_fill_gradientn(colors = palette,
                         na.value = "grey95",
                         guide = guide_colorbar(ticks = T,
                                                nbin = 50,
                                                barheight = 0.5,
                                                label = T,
                                                barwidth = 10))+
    theme_minimal()+
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 9,
        hjust = 1,
        vjust = 0.3
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.title = element_blank()
    ) +
    coord_fixed(0.5)
  
  ggsave(
    "heatmap.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 300,
    width = 40,
    height = 60
  )
  return(p)
}