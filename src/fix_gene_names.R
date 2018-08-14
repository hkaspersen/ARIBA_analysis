# Corrects the gene names found in the "cluster" column
fix_gene_names <- function(df) {
  cluster <- unique(df$cluster)
  new_names <- gsub("+", "", cluster, fixed = T)
  new_names <- gsub("_", "", new_names, fixed = T)
  new_names <- gsub("-", "", new_names, fixed = T)
  new_names <- gsub("[0-9]+", "", new_names)
  
  gene_names <- c()
  for (i in new_names) {
    p <- paste(tolower(substring(i, 1,3)), substring(i, 4), sep = "", collapse = " ")
    gene_names <- c(gene_names,p)
  }
  df2 <- data.frame(cluster,gene_names) %>%
    mutate(cluster = as.character(cluster))
  
  df <- df %>%
    left_join(df2, by = "cluster")
  
  return(df)
}