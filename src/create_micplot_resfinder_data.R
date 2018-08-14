# Create micplot data from resfinder results
# Requires data frame with mic-values called "mic"
create_micplot_resfinder_data <- function(df) {
  df <- df %>%
    select(ref, ref_name, flag, ref_ctg_change) %>%
    mutate(id = 1:n(),
           ref_ctg_change = 1,
           ref_name = gsub("(.*?).[0-9]_(.*?)$", "\\1", ref_name)) %>%
    filter(flag %in% flag_selection) %>%
    select(-flag) %>%
    spread(ref_name, ref_ctg_change) %>%
    select(-id) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate_all(funs(ifelse(. !="", ., 0))) %>%
    mutate_at(.vars = vars(-ref), funs(as.numeric(.))) %>%
    mutate(none = ifelse(rowSums(.[,-1]) > 0, 0, 1)) %>%
    group_by_at(vars(-ref)) %>%
    { mutate(ungroup(.), group = group_indices(.)) } %>%
    gather(key = "gene", value = "value", -c(ref, group)) %>%
    mutate(gene = ifelse(value == 0 & gene != "none", NA, gene),
           value = ifelse(gene == "none" & value == 0, NA, value)) %>%
    na.omit() %>%
    mutate(ref = gsub("(.*?)_resfinder", "\\1", ref)) %>%
    left_join(mic, by = "ref")
  return(df)
}