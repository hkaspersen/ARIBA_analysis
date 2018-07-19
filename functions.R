# Functions used in ARIBA data analysis

## General

### Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

### Confidence interval function
get_binCI <- function(x, n) as.numeric(setNames(binom.test(x,n)$conf.int*100,
                                                c("lwr", "upr")))

## Data import

### Identifies folder names for extracting report.tsv
folder_names <- function(filepath) {
  folders <- list.files(path = filepath, pattern = "_ariba")
  return(folders)
}

### Import ariba data from report.tsv from chosen database used in ariba
get_ariba_data <- function(filepath) {
  get_folders <- folder_names(filepath)
  
  data_list <- lapply(get_folders,
                      FUN = function(folder) {
                        read.delim(
                          paste0(filepath, "/", folder, "/", "report.tsv"),
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

### Corrects the gene names found in the "cluster" column
fix_gene_names <- function(df) {
  cluster <- unique(df$cluster)
  new_names <- gsub("+", "", cluster, fixed = T)
  new_names <- gsub("_", "", new_names, fixed = T)
  new_names <- gsub("-", "", new_names, fixed = T)
  new_names <- gsub("[0-9]+", "", new_names)
  
  gene_names <- c()
  for (i in new_names) {
    p <- paste(tolower(substring(i, 1,3)), substring(i, 4), sep = "", collapse = " ")
    gene_names <- c(gene_names,p)
  }
  df2 <- data.frame(cluster,gene_names) %>%
    mutate(cluster = as.character(cluster))
  
  df <- df %>%
    left_join(df2, by = "cluster")
  
  return(df)
}

# Vector for selection of flags in ariba report
flag_selection <- c("19","27","147","155","403",
                    "411","915","923","787","795",
                    "531","539","659","667","787","795")

### Function that returns a table with information on
### how many mutations was identified in each gene from 
### the megares database
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

### Function that returns a data frame with presence/absence of acquired genes
### from the resfinder database
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

### Function that returns a data frame with information on 
### which mutations that were found for each gene in the reports
### (only for megares data)
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

### Function that returns a table with information
### on whether or not a gene is mutated
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

### Create micplot data from resfinder results
### Requires data frame with mic-values called "mic"
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
    coord_fixed(0.5)
  
  ggsave(
    "heatmap.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 300,
    width = 40,
    height = 60
  )
  return(p)
}

### Percent acquired genes in the isolates plot
percent_acquired_genes <- function(acquired_genes_df) {
  df <- acquired_genes_df %>%
    group_by(gene, result) %>%
    count() %>%
    ungroup() %>%
    mutate(result = ifelse(result == 1, "Present", "Absent")) %>%
    spread(result, n) %>%
    mutate(total = Present + Absent,
           percent = Present/total*100) %>%
    rowwise() %>%
    mutate(lwr = get_binCI(Present, total)[1],
           upr = get_binCI(Present, total)[2])
  
  p <- ggplot(df, aes(gene, percent))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = lwr,
                      ymax = upr),
                  width = 0.4,
                  alpha = 0.4)+
    labs(x = "Genes",
         y = "Percent (%) Presence in Isolates")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    "percent_acquired_genes.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 100,
    width = 20,
    height = 20
  )
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
         dpi = 100,
         width = 20,
         height = 20,
         units = "cm")
}

### Percent plot for mutations in genes
percent_mutations <- function(mut_table_df) {
  df <- mut_table_df %>%
    group_by(gene, result) %>%
    count() %>%
    ungroup() %>%
    mutate(result = ifelse(result == 1, "Present", "Absent")) %>%
    spread(result, n, fill = 0) %>%
    mutate(total = Present + Absent,
           percent = Present/total*100) %>%
    rowwise() %>%
    mutate(lwr = get_binCI(Present, total)[1],
           upr = get_binCI(Present, total)[2],
           percent = round(percent, 1))
  
  p <- ggplot(df, aes(gene, percent))+
    geom_col(color = "black")+
    geom_text(aes(label = percent),vjust = -1)+
    geom_errorbar(aes(ymin = lwr,
                      ymax = upr),
                  width = 0.4,
                  alpha = 0.4)+
    labs(x = "Genes",
         y = "Percent (%) Presence of Mutations in Isolates")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    "percent_mutations.tiff",
    p,
    device = "tiff",
    units = "cm",
    dpi = 100,
    width = 30,
    height = 25)
}

### Presence/Absence plot for mutations
PAplot_mutations <- function(df) {
  cols <- c("1" = "#95cbee","0" = "#c9e2f6")
  
  p <- ggplot(df, aes(gene, ref, fill = result)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = cols,
                      labels = c("Mutated","Wild Type"),
                      breaks = c("1","0"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 18),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank())+
    coord_fixed(0.5)
  
  ggsave("mutation_plot.tiff",
         p,
         device = "tiff",
         dpi = 300,
         width = 40,
         height = 60,
         units = "cm")
}

### Total presence/absence plot for acquired
### genes and mutations
total_presence_absence_plot <- function(mut_table, acquired_table) {
  df <- rbind(mut_table, acquired_table) %>%
    mutate(type2 = ifelse(type == "mut", 1, 2))
  
  palette <- c("mut" = "#67a9cf", "gene" = "#ef8a62")
  
  p <- ggplot(df,
              aes(fct_reorder(gene, type2),
                  ref,
                  fill = type,
                  alpha = result)) +
    geom_tile(color = "white")+
    theme_minimal()+
    scale_fill_manual(values = palette,
                      breaks = c("mut","gene"),
                      labels = c("Mutation","Acquired Gene"))+
    labs(fill = NULL)+
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.3),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top")+
    guides(alpha = FALSE)+
    coord_fixed()
  
  ggsave("total_presence_absence_plot.tiff",
         p,
         device = "tiff",
         dpi = 100,
         width = 30,
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
