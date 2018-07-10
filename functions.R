# Functions used in ARIBA data analysis

## General

### Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

### Confidence interval function
get_binCI <- function(x, n) as.numeric(setNames(binom.test(x,n)$conf.int*100,
                                                c("lwr", "upr")))

## Data import

### Identifies folder names for extracting report.tsv
folder_names <- function(filepath,pattern) {
  folders <- list.files(path = filepath, pattern = paste0("_",pattern))
  return(folders)
}

### Import ariba data from report.tsv from chosen database used in ariba
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

## Data wrangling

### Function that returns a table with information on
### how many mutations was identified in each gene from 
### the megares database
create_no_of_mut_table <- function(df) {
  flag_selection <- c("19","27","147","155","403",
                      "411","915","923","787","795",
                      "531","539","659","667","787","795")
  df <- df %>%
    select(ref, cluster, flag, ctg_cov, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(cluster, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-c(id, flag, ctg_cov)) %>%
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

### Function that returns a data frame with presence/absence of acquired genes
### from the resfinder database
create_acquired_table <- function(df) {
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
    mutate_all(funs(ifelse(. !="", ., 0))) %>%
    gather(gene, result, -ref)
  return(df)
}

### Function that returns a data frame with information on 
### which mutations that were found for each gene in the reports
### (only for megares data)
prepare_mut_table <- function(df) {
  df <- df %>%
    select(ref, cluster, flag, ctg_cov, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    spread(cluster, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-c(id, flag, ctg_cov)) %>%
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

### Function that returns a table with information
### on whether or not a gene is mutated
create_mut_table <- function(df) {
  flag_selection <- c("19","27","147","155","403",
                      "411","915","923","787","795",
                      "531","539","659","667","787","795")
  df <- df %>%
    select(ref, cluster, flag, ctg_cov, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(cluster, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_(.*?)$", "\\1", ref)) %>%
    select(-c(id, flag, ctg_cov)) %>%
    gather(gene, mut,-ref) %>%
    group_by(gene) %>%
    filter(!any(mut == "")) %>%
    ungroup() %>%
    mutate(mut = ifelse(mut == ".", 0, 1),
           mut = as.character(mut)) %>%
    mutate(gene = gsub("[0-9]", "", gene),
           gene = gsub("-", "", gene),
           gene = gsub("_", "", gene),
           gene = gsub("+", "", gene))
  return(df)
}

### Function for selecting genes of interest and filtering the 
### columns in the data frame on the resulting vector
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

### Prepares micplot data from raw results
### Requires a data frame called "mic" with mic-values
create_micplot_megares_data <- function(df) {
  df <- df %>%
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

### Create micplot data from resfinder results
### Requires data frame with mic-values called "mic"
create_micplot_resfinder_data <- function(df) {
  df <- df %>%
    select(ref, ref_name, flag, ref_ctg_change) %>%
    mutate(id = 1:n(),
           ref_ctg_change = 1,
           ref_name = gsub("(.*?).[0-9]_(.*?)$", "\\1", ref_name)) %>%
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

## Plotting

### Heatmap on number of mutations in each gene
create_mutation_heatmap <- function(df) {
  palette <- c("#e7f0fa","#c9e2f6","#95cbee",
               "#0099dc","#4ab04a","#ffd73e",
               "#eec73a","#e29421","#f05336",
               "#ce472e")
  
  p<- ggplot(df, aes(gene, ref, fill = sum)) +
    geom_tile(color = "white") +
    labs(fill = "Mutations") +
    scale_fill_gradientn(colors = palette,
                         na.value = "grey95",
                         guide = guide_colorbar(ticks = T,
                                                nbin = 50,
                                                barheight = 0.5,
                                                label = T,
                                                barwidth = 10))+
    theme_minimal()+
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 9,
        hjust = 1,
        vjust = 0.3
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.title = element_blank()
    ) +
    coord_fixed()
  
  ggsave(
    "heamap.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 300,
    width = 30,
    height = 20
  )
  return(p)
}

### Presence/absence plot for acquired genes
PAplot_acquired_genes <- function(df) {
  cols <- c("1" = "#95cbee","0" = "#c9e2f6")
  
  p <- ggplot(df, aes(gene, ref, fill = result)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = cols,
                      labels = c("Present","Absent"),
                      breaks = c("1","0"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank())+
    coord_fixed()
  
  ggsave("acquired_genes_plot.tiff",
         p,
         device = "tiff",
         dpi = 300,
         width = 20,
         height = 20,
         units = "cm")
}

### Presence/Absence plot for mutations
PAplot_mutations <- function(df) {
  cols <- c("1" = "#95cbee","0" = "#c9e2f6")
  
  p <- ggplot(df, aes(gene, ref, fill = mut)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = cols,
                      labels = c("Mutated","Wild Type"),
                      breaks = c("1","0"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank())+
    coord_fixed()
  
  ggsave("mutation_plot.tiff",
         p,
         device = "tiff",
         dpi = 300,
         width = 20,
         height = 20,
         units = "cm")
}

### Create micplot on selected antimicrobial
### data from "create_micplot_data"
create_micplot <- function(df, am, plotname) {
  factor_breaks <- c(0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 
                     128, 256, 512, 1024, 2048)
  
  factor_levels <- c("0.015", "0.03", "0.06", "0.12", "0.25", "0.5", "1", "2", "4", 
                     "8", "16", "32", "64", "128", "256", "512", "1024", "2048")
  
  p1 <- ggplot(df, aes(factor(group),gene,fill = factor(group)))+
    geom_point(pch = 21, size = 3)+
    guides(fill = FALSE)+
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = 7))
  
  
  p2 <- ggplot(df, aes(factor(group),am, fill = factor(group)))+
    geom_boxplot()+
    scale_y_continuous(labels = factor_levels,
                       breaks = factor_breaks,
                       trans = "log10")+
    theme(axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank())+
    guides(fill = FALSE)
  
  p3 <- plot_grid(
    p2,
    p1,
    align = "v",
    nrow = 2,
    rel_heights = c(3 / 5, 2 / 5)
  )
  
  ggsave(
    paste0(plotname, ".tiff"),
    p3,
    device = "tiff",
    dpi = 300,
    width = 20,
    height = 30,
    units = "cm"
  )
}
