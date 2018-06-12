library(tidyverse)

## Functions

# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

# Identifies folder names for extracting report.tsv
folder_names <- function(filepath,pattern) {
  folders <- list.files(path = filepath, pattern = paste0("_",pattern))
  return(folders)
}

# Import ariba data from report.tsv from chosen database used in ariba
get_ariba_data <- function(path, database) {
  get_folders <- folder_names(path, database)
  
  data_list <- lapply(get_folders,
                      FUN = function(folder) {
                        read.delim(
                          paste0(database, "/", folder, "/", "report.tsv"),
                          stringsAsFactors = F,
                          header = TRUE,
                          sep = "\t"
                        )
                      })
  
  names(data_list) <- get_folders
  data <- bind_rows(data_list, .id = "ref")
  return(data)
}

# Summarizes data from report.tsv
summarise_data <- function(df) {
  df <- df %>%
    select(ref, cluster, flag, ctg_cov, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    spread(cluster, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-id)
  df[df == ""] <- "0"
  df[df == "."] <- "1"
  return(df)
}

# Split all mutations in columns to single 1/0 columns
AddCol <- function(df, col_name) {
  # split rows by delimeters
  string_to_proc <- df %>% select(!!col_name) %>%
    unlist() %>% str_split(regex("\\, |\\,")) 
  # find unique entries
  unique_strings <- string_to_proc %>%
    unlist() %>% unique()
  # construct names of the new columns
  cols_names <- paste(col_name, unique_strings, sep = "_")
  # construct 0/1-content columns for each unique entry
  cols_content <- sapply(function(i) {
    as.integer(unlist(lapply(function(Z) any(Z %in% unique_strings[i]), 
                             X = string_to_proc)))
  }, X = seq_along(unique_strings))
  res <- data.frame(cols_content)
  names(res) <- cols_names
  return(res)
}

# Create data frame report from megares data
create_mut_report <- function(df) {
  flag_selection <- c("19","27","147","155","403",
                      "411","915","923","787","795",
                      "531","539","659","667","787","795")
  
  df <- df %>%
    select(ref, cluster, flag, ref_ctg_change) %>%
    filter(ref_ctg_change != ".",
           flag %in% flag_selection) %>%
    mutate(id = 1:n()) %>%
    spread(cluster, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, flag))
  
  col_names <- colnames(df)
  col_names <- col_names[-1]
  
  cols_to_add <- sapply(function(i) {AddCol(df = df, col_name = col_names[i])}, 
                        X = seq_along(col_names)) %>% 
    bind_cols()
  
  df <- df %>%
    bind_cols(cols_to_add)
  
  df_complete <- df[,!names(df) %in% col_names]
  
  return(df_complete)
}

# Creates data frame report from resfinder data
create_AG_report <- function(df) {
  flag_selection <- c("19","27","147","155","403",
                      "411","915","923","787","795",
                      "531","539","659","667","787","795")
  
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
    mutate_all(funs(ifelse(. !="", ., 0)))
  return(df)
}

## Import files and create reports

# Megares
mut_data <- get_ariba_data("D:\\R-Projects\\Ariba_analysis\\megares", "megares")
mut_report <- create_mut_report(mut_data)

# Resfinder
gene_data <- get_ariba_data("D:\\R-Projects\\Ariba_analysis\\resfinder", "resfinder")
gene_report <- create_AG_report(gene_data)

# Statistics

test <- mut_report %>%
  gather(gene_mut, value, -ref) %>%
  mutate(gene = gsub("^(.*?)_(.*?)$", "\\1", gene_mut),
         gene = gsub("[0-9]", "", gene),
         gene = gsub("^(.*?)_(.*?)$", "\\1", gene),
         gene = gsub("^(.*?)_$", "\\1", gene),
         gene = gsub("+", "", gene),
         gene = gsub("-", "", gene))



stat_df <- data.frame(t(mut_report))
colnames(stat_df) <- mut_report$ref
stat_df$ref <- row.names(stat_df)
stat_df <- stat_df[-1,]
row.names(stat_df) <- 1:length(stat_df[,1])
stat_df_complete <- stat_df[,c(11,1:10)]

stat_report <- stat_df_complete %>%
  mutate_at(.vars = vars(-ref),
            funs(as.numeric(as.character(.)))) %>%
  mutate(total = ncol(.)-1,
         no_samples_mut = rowSums(.[,2:11]))
  









