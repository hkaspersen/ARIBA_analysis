# ARIBA analysis script

## Information
This is an R script designed to work in a shell bash environment. The 
script produces multiple plots and files with information related to 
mutations in intrinsic genes and the presence/absence of acquired genes.

The script is designed to work with files from the Bifrost pipeline, 
found [here](https://github.com/NorwegianVeterinaryInstitute/Bifrost), 
run on the MegaRes and ResFinder databases, in separate folders.

## Usage
Make sure that your .Rprofile file have the correct library path before 
beginning. The .Rprofile file is located in your home folder, and you 
get there by typing **cd**. Type in **nano .Rprofile**, and paste in the 
following:

```
path <- "/work/projects/nn9305k/lib/R/"

.libPaths(c(path, .libPaths()))
```

Load up R in your terminal: 
```
module load R/3.5.0
```
Then run the following command:
```
Rscript /work/projects/nn9305k/vi_src/ARIBA_analysis/ariba_analysis.R 
megares_reports_path 
resfinder_reports_path 
output_dir_path acquired_genes_file intrinsic_genes_file 
heatmap_selection
```
- **megares_reports_path**: Full path to the folder holding the megares 
reports
- **resfinder_reports_path**: Full path to the folder holding the 
resfinder reports
- **output_dir_path**: Full path to the folder where output files will 
be placed
- **acquired_genes_file**: Path to file holding the list of acquired 
genes of interest. Note: It is enough to write "qnr" only to get all qnr 
genes. Same rule applies to all genes. To get all genes reported by 
resfinder, simply input "all" in this file
- **intrinsic_genes_file**: Path to file holding the list of intrisic 
genes of interest. Same rules as above.
- **heatmap_selection**: TRUE/FALSE, whether or not the heatmap should 
be generated. Note that this will fail if there are no Qnr genes present 
in the dataset (to be fixed).

## File descriptions
Plots

- **acquired_stats.svg**: A figure showing basic statistics for the 
chosen acquired genes, including the percent of isolates with these genes 
present, and a 95 % confidence interval.
- **mut_stats.svg**: Similar to acquired_stats.svg, but only for the 
chosen intrinsic genes on interest. The figure shows the percent of isolates 
with mutations in the intrinsic genes. The script includes mutation 
position correction for *gyrA*, *gyrB*, *parC* and *parE*.
- **heatmap.svg**: A heatmap over mutations in *gyrA*, *gyrB*, *parC* 
and *parE*, along with presence/absence of qnr-genes. 

Tables

- **acquired_report**: Contains present/absent data on selected acquired 
genes for each isolate included. One row per isolate.
- **mut_report**: Similar to acquired_report, except with mutations in 
intrinsic genes.
- **mut_flags**: Contains a list over all mutations reported by ARIBA 
MegaRes, and their corresponding flag. The "flag_result" column reports 
whether or not the flag passed the filtering (1 = passed).
- **acquired_flags**: Same as above, but for acquired genes.
- **mut_quant**: Contains information about which mutation was found in 
*gyrA*, *gyrB*, *parC*, and *parE*, and also how many mutations were 
found in each respective gene (mut_gyrA etc). Each row is an isolate.
- **no_of_mut**: Contains the number of isolates for each number of 
mutations group. Each row corresponds to one group, where "n" is the 
amount of isolates with the corresponding number of mutations listed in 
the other columns.
- **mut_combinations**: Contains the number of isolates for each 
mutation in the genes *gyrA*, *gyrB*, *parC*, and *parE*. The actual 
mutations are listed here.
