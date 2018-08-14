# Total presence/absence plot for acquired genes and mutations
total_presence_absence_plot <- function(mut_table, acquired_table) {
  df <- rbind(mut_table, acquired_table) %>%
    mutate(type2 = ifelse(type == "mut", 1, 2))
  
  palette <- c("mut" = "#67a9cf", "gene" = "#ef8a62")
  
  p <- ggplot(df,
              aes(fct_reorder(gene, type2),
                  ref,
                  fill = type,
                  alpha = result)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = palette,
                      breaks = c("mut","gene"),
                      labels = c("Mutation","Acquired Gene"))+
    labs(fill = NULL)+
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.3),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top")+
    guides(alpha = FALSE)+
    coord_fixed()
  
  ggsave("total_presence_absence_plot.tiff",
         p,
         device = "tiff",
         dpi = 100,
         width = 30,
         height = 20,
         units = "cm")
}