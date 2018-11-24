#!/usr/bin/env Rscript

####################################################
################## ARIBA ANALYSIS ##################
####################################################

# This script takes multiple input from the software
# ARIBA (github.com/sanger-pathogens/ariba) and
# create summary reports and visualizations. One
# can select specific analyses based on the data you
# have: AMR, virulence, MLST or plasmids.

# Author: Håkon Kaspersen

####################################################
####################################################
####################################################

# ------------------- Libraries --------------------
suppressPackageStartupMessages(pacman::p_load(optparse))
parser <- OptionParser(usage = "Usage: %prog [options]")

# Create command line options
parser <- add_option(parser,
                     c("-a", "--amr"),
                     action = "store_true",
                     help = "Run AMR gene analysis.")
parser <- add_option(parser,
                     c("-u", "--mut"),
                     action = "store",
                     help = "Location of intrinsic gene reports.")
parser <- add_option(parser,
                     c("-q", "--acq"),
                     action = "store",
                     help = "Location of acquired gene reports.")
parser <- add_option(parser,
                     c("-i", "--intrinsic"),
                     action = "store",
                     help = "List of intrinsic genes of interest. Type 'all' for including all reported genes.")
parser <- add_option(parser,
                     c("-c", "--acquired"),
                     action = "store",
                     help = "List of acquired genes of interest. Type 'all' for including all reported genes.")
parser <- add_option(parser,
                     c("-v", "--vir"),
                     action = "store",
                     help = "Location of ARIBA virulence reports.")
parser <- add_option(parser,
                     c("-m", "--mlst"),
                     action = "store",
                     help = "Location of ARIBA MLST reports.")
parser <- add_option(parser,
                     c("-p", "--plasmid"),
                     action = "store",
                     help = "Location of ARIBA plasmid reports.")
parser <- add_option(parser,
                     c("-o", "--output"),
                     action = "store",
                     help = "Output directory location. One folder for each analysis will be created at given location.")
opt <- parse_args(parser)

# Check if output folder is specified
if (is.null(opt$output)) {
  print("Please specify an output directory.")
  stop()
}

## ------------------- Tracks ----------------------
## AMR gene track
if (!is.null(opt$amr)) {
  print(paste0(
    "Running AMR gene summary analysis. Intrinsic gene reports location: ",
    opt$u,
    ". Acquired gene reports location: ",
    opt$q,
    ". Output location: ",
    opt$o))
  system(paste("Rscript amr_script.R",
               opt$u, 
               opt$q, 
               opt$o, 
               opt$c,
               opt$i))
}

## Virulence gene track
if (!is.null(opt$vir)) {
  print(paste0(
    "Running virulence gene summary analysis. Reports location: ",
    opt$vir,
    ". Output location: ",
    opt$o))
  system(paste("Rscript vir_script.R", opt$v, opt$o))
}

## MLST track
if (!is.null(opt$mlst)) {
  print(paste0(
    "Running MLST summary analysis. Reports location: ", 
    opt$mlst,
    ". Output location: ",
    opt$o))
  system(paste("Rscript mlst_script.R", opt$m, opt$o))
}

## Plasmid typing track
if (!is.null(opt$plasmid)) {
  print(paste0(
    "Running plasmid summary analysis. Reports location: ",
    opt$plasmid,
    ". Output location: ",
    opt$o))
  system(paste("Rscript plasmid_script.R", opt$p, opt$o))
}
