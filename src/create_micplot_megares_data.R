# Prepares micplot data from raw results
# Requires a data frame called "mic" with mic-values
create_micplot_megares_data <- function(df) {
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
    mutate(mut = ifelse(mut == ".", 0, 1)) %>%
    spread(gene, mut) %>%
    mutate(none = ifelse(rowSums(.[,-1]) > 0, 0, 1)) %>%
    group_by_at(vars(-ref)) %>%
    { mutate(ungroup(.), group = group_indices(.)) } %>%
    gather(key = "gene", value = "value", -c(ref, group)) %>%
    mutate(gene = ifelse(value == 0 & gene != "none", NA, gene),
           value = ifelse(gene == "none" & value == 0, NA, value)) %>%
    na.omit() %>%
    left_join(mic, by = "ref")
  return(df)
}