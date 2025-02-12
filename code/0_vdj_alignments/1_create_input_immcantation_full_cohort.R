#!/usr/bin/env Rscript

# START ----
# Workflow for processing combined analysis object
#  - Prior to running this code set working directory to 'impact-analysis/code/0_vdj_alignments'
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings)) # for FASTAs
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(impactSingleCellToolkit))


# INPUT and OUTPUT Directories
param_file_fh = "../../input/input_full_cohort_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads

# OUTPUT Results Directories
vdj_bcr_dir <- "../../../Results/vdj_bcr_full_cohort"
vdj_tcr_dir <- "../../../Results/vdj_tcr_full_cohort"
results_dir <- "../../../Results/input_immcantation_full_cohort"
if(!file.exists(vdj_bcr_dir)){dir.create(vdj_bcr_dir,recursive = TRUE)}
if(!file.exists(vdj_tcr_dir)){dir.create(vdj_tcr_dir,recursive = TRUE)}
if(!file.exists(results_dir)){dir.create(results_dir,recursive = TRUE)}

# Separate TCR and BCR metadata
metadata_bcr <- metadata[metadata$type %in% "vdj_b",]
metadata_tcr <- metadata[metadata$type %in% "vdj_t",]

# Create symbolic links of select sample results to VDJ sub-directories

## BCR subset
bcr_list        <- metadata_bcr$results_directory_path
names(bcr_list) <- metadata_bcr$sample
bcr_list        <- as.list(bcr_list)

for (sample in names(bcr_list)) {
  sample_dir  <- bcr_list[[sample]]
  symlink_dir <- file.path(vdj_bcr_dir,sample)
  
  # Confirm filtered annotations file is not empty
  filtered_results_exist <- TRUE
  files        <- list.files(sample_dir)
  filtered_csv <- files[grep("\\.csv",files)]
  filtered_csv <- filtered_csv[grep("filtered",filtered_csv)]
  info         <- file.info(file.path(sample_dir,filtered_csv))
  if(info$size == as.integer(0)){
    filtered_results_exist <- FALSE
  }
  
  # Create the symbolic link
  if(!file.exists(symlink_dir)){
    if(filtered_results_exist){
      system(paste("ln -s", sample_dir, symlink_dir))
    }
  }
}

## TCR subset
tcr_list        <- metadata_tcr$results_directory_path
names(tcr_list) <- metadata_tcr$sample
tcr_list        <- as.list(tcr_list)

for (sample in names(tcr_list)) {
  sample_dir  <- tcr_list[[sample]]
  symlink_dir <- file.path(vdj_tcr_dir,sample)
  
  # Confirm filtered annotations file is not empty
  filtered_results_exist <- TRUE
  files        <- list.files(sample_dir)
  filtered_csv <- files[grep("\\.csv",files)]
  filtered_csv <- filtered_csv[grep("filtered",filtered_csv)]
  info         <- file.info(file.path(sample_dir,filtered_csv))
  if(info$size == as.integer(0)){
    filtered_results_exist <- FALSE
  }
  
  # Create the symbolic link
  if(!file.exists(symlink_dir)){
    if(filtered_results_exist){
      system(paste("ln -s", sample_dir, symlink_dir))
    }
  }
}


# Generate input to immcantation

## BCR 
generateInputFilesImmcantation(vdj.results.dir = vdj_bcr_dir,
                               contig.annotation.file = "filtered_contig_annotations.csv",
                               contig.fasta.file = "filtered_contig.fasta",
                               output.file.prefix = "input_bcr",
                               output.dir = results_dir)

## TCR 
generateInputFilesImmcantation(vdj.results.dir = vdj_tcr_dir,
                               contig.annotation.file = "filtered_contig_annotations.csv",
                               contig.fasta.file = "filtered_contig.fasta",
                               output.file.prefix = "input_tcr",
                               output.dir = results_dir)
