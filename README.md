# ARIBA_analysis

## Information
This is an R script designed to work in a shell bash environment. The 
script produces multiple plots and files with information related to 
mutations in intrinsic genes and the presence/absence of acquired genes.

## Usage
First, load up R in your terminal: **module load R/3.5.0**

Then run the following command:

**Rscript ariba_analysis.R megares_reports_path resfinder_reports_path 
output_dir_path acquired_genes_file intrinsic_genes_file**

- megares_reports_path: Full path to the folder holding the megares 
reports
- resfinder_reports_path: Full path to the folder holding the resfinder 
reports
- output_dir_path: Full path to the folder where output files will be 
placed
- acquired_genes_file: Path to file holding the list of acquired genes 
of interest. Note: It is enough to write "qnr" only to get all qnr 
genes. Same rule applies to all genes. To get all genes reported by 
resfinder, simply input "all" in this file
- intrinsic_genes_file: Path to file holding the list of intrisic genes 
of interest. Same rules as above.
