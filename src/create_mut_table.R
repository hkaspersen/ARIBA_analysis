# Function that returns a table with information
# on whether or not a gene is mutated
create_mut_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    group_by(gene) %>%
    filter(!any(result == "")) %>%
    ungroup() %>%
    mutate(result = ifelse(result == ".", 0, 1),
           result = as.character(result)) %>%
    mutate(gene = gsub("[0-9]", "", gene),
           gene = gsub("-", "", gene),
           gene = gsub("_", "", gene),
           gene = gsub("+", "", gene),
           type = "mut")
  return(df)
}