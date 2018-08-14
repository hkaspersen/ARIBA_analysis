# Barplot of percent acquired genes
percent_acquired_genes <- function(acquired_genes_df) {
  df <- acquired_genes_df %>%
    group_by(gene, result) %>%
    count() %>%
    ungroup() %>%
    mutate(result = ifelse(result == 1, "Present", "Absent")) %>%
    spread(result, n) %>%
    mutate(total = Present + Absent,
           percent = Present/total*100) %>%
    rowwise() %>%
    mutate(lwr = get_binCI(Present, total)[1],
           upr = get_binCI(Present, total)[2])
  
  p <- ggplot(df, aes(gene, percent))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = lwr,
                      ymax = upr),
                  width = 0.4,
                  alpha = 0.4)+
    labs(x = "Genes",
         y = "Percent (%) Presence in Isolates")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    "percent_acquired_genes.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 100,
    width = 20,
    height = 20
  )
}