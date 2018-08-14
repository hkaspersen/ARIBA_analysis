# Function that returns a table with information on
# how many mutations was identified in each gene from 
# the megares database
create_no_of_mut_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, mut,-ref) %>%
    group_by(gene) %>%
    filter(!any(mut == "")) %>%
    ungroup() %>%
    mutate(mut = ifelse(mut == ".", NA, mut),
           entries = sapply(strsplit(.$mut, ","), FUN = function(x) {length(x)})) %>%
    separate(mut, into = as.character(c(1:max(.$entries)))) %>%
    gather(no_mut, mut,-ref,-gene) %>%
    mutate(mut = ifelse(mut == "", NA, ifelse(!is.na(mut), 1, mut)),
           mut = as.numeric(mut)) %>%
    spread(no_mut, mut) %>%
    mutate(sum = rowSums(subset(., select = -c(ref, gene)), na.rm = TRUE)) %>%
    select(ref, gene, sum) %>%
    mutate(gene = gsub("[0-9]", "", gene),
           gene = gsub("-", "", gene),
           gene = gsub("_", "", gene),
           gene = gsub("+", "", gene))
  return(df)
}