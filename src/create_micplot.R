# Creates micplot
create_micplot <- function(df, am, plotname) {
  factor_breaks <- c(0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 
                     128, 256, 512, 1024, 2048)
  
  factor_levels <- c("0.015", "0.03", "0.06", "0.12", "0.25", "0.5", "1", "2", "4", 
                     "8", "16", "32", "64", "128", "256", "512", "1024", "2048")
  
  p1 <- ggplot(df, aes(factor(group),gene,fill = factor(group)))+
    geom_point(pch = 21, size = 3)+
    guides(fill = FALSE)+
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = 7))
  
  
  p2 <- ggplot(df, aes(factor(group),am, fill = factor(group)))+
    geom_boxplot()+
    scale_y_continuous(labels = factor_levels,
                       breaks = factor_breaks,
                       trans = "log10")+
    theme(axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank())+
    guides(fill = FALSE)
  
  p3 <- plot_grid(
    p2,
    p1,
    align = "v",
    nrow = 2,
    rel_heights = c(3 / 5, 2 / 5)
  )
  
  ggsave(
    paste0(plotname, ".tiff"),
    p3,
    device = "tiff",
    dpi = 300,
    width = 20,
    height = 30,
    units = "cm"
  )
}
