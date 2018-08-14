# Percent plot for mutations in genes
percent_mutations <- function(mut_table_df) {
  df <- mut_table_df %>%
    group_by(gene, result) %>%
    count() %>%
    ungroup() %>%
    mutate(result = ifelse(result == 1, "Present", "Absent")) %>%
    spread(result, n, fill = 0) %>%
    mutate(total = Present + Absent,
           percent = Present/total*100) %>%
    rowwise() %>%
    mutate(lwr = get_binCI(Present, total)[1],
           upr = get_binCI(Present, total)[2],
           percent = round(percent, 1))
  
  p <- ggplot(df, aes(gene, percent))+
    geom_col(color = "black")+
    geom_text(aes(label = percent),vjust = -1)+
    geom_errorbar(aes(ymin = lwr,
                      ymax = upr),
                  width = 0.4,
                  alpha = 0.4)+
    labs(x = "Genes",
         y = "Percent (%) Presence of Mutations in Isolates")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    "percent_mutations.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 100,
    width = 30,
    height = 25)
}