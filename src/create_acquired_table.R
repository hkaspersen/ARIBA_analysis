# Function that returns a data frame with presence/absence of acquired genes
# from the resfinder database
create_acquired_table <- function(df) {
  df <- df %>%
    select(ref, ref_name, flag, ref_ctg_change) %>%
    mutate(id = 1:n(),
           ref_ctg_change = 1,
           ref_name = gsub("(.*?).[0-9]_(.*?)$", "\\1", ref_name),
           ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    filter(flag %in% flag_selection) %>%
    select(-flag) %>%
    spread(ref_name, ref_ctg_change) %>%
    select(-id) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate_all(funs(ifelse(. !="", ., 0))) %>%
    gather(gene, result, -ref) %>%
    mutate(type = "gene")
  return(df)
}