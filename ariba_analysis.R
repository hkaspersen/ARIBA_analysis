# Running script for ARIBA analysis

# ------------------------------------ Libraries ----------------------------------------

library(tidyverse)
library(cowplot)
library(gridExtra)
library(grid)
library(forcats)

# -------------------------------- Import functions -------------------------------------

source("functions.R")

# ---------------------------------- Import data ----------------------------------------
mut_data <- get_ariba_data("D:\\R-Projects\\Ariba_analysis\\ariba_run_megares")
gene_data <- get_ariba_data("D:\\R-Projects\\Ariba_analysis\\resfinder", "resfinder")

# ------------------------------- Run data wrangling ------------------------------------

### Fix gene names in data frame
clean_mut_data <- fix_gene_names(mut_data)

### Amount of mutations per gene
quantified_mutations <- create_no_of_mut_table(clean_mut_data)

### Presence/Absence of mutations
mut_table <- create_mut_table(clean_mut_data)

### Information on which mutations
actual_mutations_all <- prepare_mut_table(clean_mut_data)

actual_mutations_filtered <- select_genes(actual_mutations_all,
                                          c("gyra","gyrb","parc","pare",
                                            "marr","mara","tolc","ompf"))
### Acquired genes
acquired_table <- create_acquired_table(gene_data)

### Prepare micplot data
mic <- read.table("mic_values.txt", sep = "\t", stringsAsFactors = F, header = T)

micplot_megares_data <- create_micplot_megares_data(clean_mut_data)
micplot_resfinder_data <- create_micplot_resfinder_data(gene_data)

# ---------------------------------- Plotting -------------------------------------------

### Quantified mutations heatmap
create_mutation_heatmap(quantified_mutations)

### Presence/Absence of mutations figure
PAplot_mutations(mut_table)

### Acquired genes plot
PAplot_acquired_genes(acquired_table)

### Percent plots
percent_mutations(mut_table)
percent_acquired_genes(acquired_table)

### Total presence/absence plot
total_presence_absence_plot(mut_table, acquired_table)

### Micplots

#### Megares
create_micplot(micplot_megares_data, micplot_megares_data$T_CIP, "micplot_cip_MR")
create_micplot(micplot_megares_data, micplot_megares_data$U_NAL, "micplot_nal_MR")

#### Resfinder
create_micplot(micplot_resfinder_data, micplot_resfinder_data$T_CIP, "micplot_cip_RF")
create_micplot(micplot_resfinder_data, micplot_resfinder_data$U_NAL, "micplot_nal_RF")
