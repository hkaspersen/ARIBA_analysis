library(tidyverse)

## Functions

# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

# Confidence interval function
get_binCI <- function(x, n) as.numeric(setNames(binom.test(x,n)$conf.int*100,
                                                c("lwr", "upr")))

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
  

gene_stat_df <- data.frame(t(gene_report))
colnames(gene_stat_df) <- gene_report$ref
gene_stat_df$ref <- row.names(gene_stat_df)
gene_stat_df <- gene_stat_df[-1,]
row.names(gene_stat_df) <- 1:length(gene_stat_df[,1])
gene_stat_df_complete <- gene_stat_df[,c(11,1:10)]

gene_report_df <- gene_stat_df_complete %>%
  mutate_at(.vars = vars(-ref),
            funs(as.numeric(as.character(.)))) %>%
  mutate(total = ncol(.)-1,
         no_samples = rowSums(.[,2:11]),
         percent_genes = no_samples/total *100) %>%
  rowwise() %>%
  mutate(lwr = get_binCI(no_samples, total)[1],
         upr = get_binCI(no_samples, total)[2])


ggplot(gene_report_df, aes(ref, percent_genes)) +
  geom_col() +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                alpha = 0.3,
                width = 0.4)+
  labs(y = "Percent of Samples")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

# -------------------------------------------------------------------------------------
## Testing for presence of mutations in genes only

heatmap_df <- mut_data %>%
  select(ref, cluster, flag, ctg_cov, ref_ctg_change) %>%
  mutate(id = 1:n()) %>%
  spread(cluster, ref_ctg_change) %>%
  group_by(ref) %>%
  summarise_all(funs(func_paste)) %>%
  mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
  select(-c(id, flag, ctg_cov)) %>%
  gather(gene, mut,-ref) %>%
  group_by(gene) %>%
  filter(!any(mut == "")) %>%
  ungroup() %>%
  mutate(mut = ifelse(mut == ".", NA, mut)) %>%
  separate(mut, into = as.character(c(1:35))) %>%
  gather(no_mut, mut,-ref,-gene) %>%
  mutate(mut = ifelse(mut == "", NA, ifelse(!is.na(mut), 1, mut)),
         mut = as.numeric(mut)) %>%
  spread(no_mut, mut) %>%
  mutate(sum = rowSums(subset(., select = -c(ref, gene)), na.rm = TRUE)) %>%
  select(ref, gene, sum)

ggplot(heatmap_df, aes(gene, ref, fill = sum)) +
  geom_tile(color = "white") +
  labs(fill = "Mutations") +
  scale_fill_distiller(palette = 1)+
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 9,
      hjust = 1,
      vjust = 0.3
    ),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_fixed()

ggsave(
  "heamap.tiff",
  p,
  device = "tiff",
  units = "cm",
  dpi = 300,
  width = 40,
  height = 30
)

?scale_fill_distiller






