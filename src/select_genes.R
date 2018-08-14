# Function for selecting genes of interest and filtering the 
# columns in the data frame on the resulting vector
select_genes <- function(df, gene_string) {
  names_df <- names(df)
  x <- "ref"
  for (i in names_df) {
    for (j in gene_string) {
      if (str_detect(i, regex(j, ignore_case = T)) == TRUE) {
        x <- c(x, i)
      }
    }
  }
  df <- df %>%
    select(one_of(x))
  return(df)
}