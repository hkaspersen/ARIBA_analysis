# Running script for ARIBA analysis

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

micplot_mut_data <- create_micplot_data(mut_data)

# ---------------------------------- Plotting -------------------------------------------

### Quantified mutations heatmap
quant_heat <- create_mutation_heatmap(quantified_mutations)

ggsave("quantified_mutations_heatmap.tiff",
       quant_heat,
       device = "tiff",
       dpi = 300,
       width = 35,
       height = 25,
       units = "cm")

### Acquired genes plot

acquired_plot <- PAplot_acquired_genes(acquired_genes)

ggsave("acquired_genes_plot.tiff",
       acquired_plot,
       device = "tiff",
       dpi = 300,
       width = 20,
       height = 20,
       units = "cm")

### Micplots

cip_plot <- create_micplot(micplot_data, micplot_data$T_CIP)
nal_plot <- create_micplot(micplot_data, micplot_data$U_NAL)

ggsave("cip_plot.tiff",
       cip_plot,
       device = "tiff",
       dpi = 300,
       width = 20,
       height = 30,
       units = "cm")

ggsave("nal_plot.tiff",
       nal_plot,
       device = "tiff",
       dpi = 300,
       width = 20,
       height = 30,
       units = "cm")



