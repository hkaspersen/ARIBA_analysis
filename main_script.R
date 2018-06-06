library(tidyverse)

## Functions

# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

func_paste2 <- function(x) paste(unique(sum(x), collapse = ","))

# Identifies folder names for extracting report.tsv
folder_names <- function(filepath,pattern) {
  folders <- list.files(path = filepath, pattern = paste0("_",pattern))
  return(folders)
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


## Import files

# Megares

mut_folders <- folder_names("D:\\R-Projects\\Ariba_analysis\\megares", "megares")

mut_data_list <- lapply(mut_folders,
       FUN = function(folder) {
         read.delim(
           paste0("megares/", folder, "/", "report.tsv"),
           stringsAsFactors = F,
           header = TRUE,
           sep = "\t"
         )
       })

names(mut_data_list) <- mut_folders

mut_data <- bind_rows(mut_data_list, .id = "ref")

# Resfinder

acquired_folders <- folder_names("D:\\R-Projects\\Ariba_analysis\\resfinder", "resfinder")

acquired_data_list <- lapply(acquired_folders,
                        FUN = function(folder) {
                          read.delim(
                            paste0("resfinder/", folder, "/", "report.tsv"),
                            stringsAsFactors = F,
                            header = TRUE,
                            sep = "\t"
                          )
                        })

names(acquired_data_list) <- acquired_folders

acquired_data <- bind_rows(acquired_data_list, .id = "ref")

## Data wrangling

# Megares

flag_selection <- c("19","27","147","155","403",
                    "411","915","923","787","795",
                    "531","539","659","667","787","795")

test <- mut_data %>%
  select(ref, cluster, flag, ref_ctg_change) %>%
  filter(ref_ctg_change != ".",
         flag %in% flag_selection) %>%
  mutate(id = 1:n()) %>%
  spread(cluster, ref_ctg_change) %>%
  group_by(ref) %>%
  summarise_all(funs(func_paste)) %>%
  select(-c(id, flag))

col_names <- colnames(test)
col_names <- col_names[-1]

cols_to_add <- sapply(function(i) {AddCol(df = test, col_name = col_names[i])}, 
                      X = seq_along(col_names)) %>% 
  bind_cols()

test <- test %>%
  bind_cols(cols_to_add)

test_complete <- test[,!names(test) %in% col_names]

mut_report <- test_complete %>%
  gather(key, value, -ref)

# Resfinder

acquired_data_test <- acquired_data %>%
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













