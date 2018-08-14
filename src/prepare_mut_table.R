# Function that returns a data frame with information on 
# which mutations that were found for each gene in the reports
# (only for megares data)
prepare_mut_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    filter(flag %in% flag_selection) %>%
    mutate(id = 1:n()) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, mut,-ref) %>%
    mutate(gene = gsub("[0-9]", "", gene),
           gene = gsub("-", "", gene),
           gene = gsub("_", "", gene),
           gene = gsub("+", "", gene)) %>%
    group_by(gene) %>%
    filter(!any(mut == "")) %>%
    ungroup() %>%
    mutate(mut = ifelse(mut == ".", NA, mut),
           id = 1:n()) %>%
    spread(gene, mut) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    select(-id)
  return(df)
}