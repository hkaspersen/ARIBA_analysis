# Running script for ARIBA analysis

# ------------------------------------ Libraries ----------------------------------------

library(tidyverse)
library(cowplot)
library(gridExtra)
library(grid)

# -------------------------------- Import functions -------------------------------------

source("functions.R")

# ---------------------------------- Import data ----------------------------------------

mut_data <- get_ariba_data("D:\\R-Projects\\Ariba_analysis\\megares", "megares")
gene_data <- get_ariba_data("D:\\R-Projects\\Ariba_analysis\\resfinder", "resfinder")

# ------------------------------- Run data wrangling ------------------------------------

### Amount of mutations per gene
quantified_mutations <- create_no_of_mut_table(mut_data)

### Information on which mutations
actual_mutations_all <- prepare_mut_table(mut_data)

actual_mutations_filtered <- select_genes(actual_mutations_all,
                                          c("gyra","gyrb","parc","pare",
                                            "marr","mara","tolc","ompf"))
### Acquired genes
acquired_genes <- create_acquired_table(gene_data)

### Prepare micplot data
mic <- read.table("mic_values.txt", sep = "\t", stringsAsFactors = F, header = T)

micplot_megares_data <- create_micplot_megares_data(mut_data)
micplot_resfinder_data <- create_micplot_resfinder_data(gene_data)

# ---------------------------------- Plotting -------------------------------------------

### Quantified mutations heatmap
quant_heat <- create_mutation_heatmap(quantified_mutations)

### Acquired genes plot

acquired_plot <- PAplot_acquired_genes(acquired_genes)

### Micplots

#### Megares

cip_plot_megares <- create_micplot(micplot_megares_data, micplot_megares_data$T_CIP, "micplot_cip_MR")
nal_plot_megares <- create_micplot(micplot_megares_data, micplot_megares_data$U_NAL, "micplot_nal_MR")

#### Resfinder

cip_plot_resfinder <- create_micplot(micplot_resfinder_data, micplot_resfinder_data$T_CIP, "micplot_cip_RF")
nal_plot_resfinder <- create_micplot(micplot_resfinder_data, micplot_resfinder_data$U_NAL, "micplot_nal_RF")
