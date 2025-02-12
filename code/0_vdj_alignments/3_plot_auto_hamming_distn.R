#!/usr/bin/env Rscript

# START ----
# Code to plot hamming distance disctribution across immcantation results
#  - Prior to running this code set working directory to 'impact-analysis/code/0_vdj_alignments'

# Immcantation
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dowser))
suppressPackageStartupMessages(library(airr))

suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(impactSingleCellToolkit))

# If running in the terminal add a setwd line below to the folder 'impact-analysis/code/0_vdj_alignments'
#setwd("")

# INPUT and OUTPUT Directories
param_file_fh = "../../input/input_full_cohort_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
THREADS       = params$threads

# OUTPUT Results Directories
plot_dir     = params$plot_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}


# Load input results
bcr_heavy_fh <- "../../../Results/results_immcantation_full_cohort/BCR_CSFPB/BCR_CSFPB_heavy_germ-pass.tsv"
bcr_light_fh <- "../../../Results/results_immcantation_full_cohort/BCR_CSFPB/BCR_CSFPB_light_germ-pass.tsv"

db.heavy <- read_rearrangement(bcr_heavy_fh) %>% as.data.frame()
db.light <- read_rearrangement(bcr_light_fh) %>% as.data.frame()

# Merge chains
bcr_df <- mergeVDJresults(db.heavy,db.light,umi.thresh = 3,assay = "bcr")

# format VDJ unique cell IDs to match scRNA-Seq
formatUniqueIDcolum <- function(uni_id_column,sep = ":") {
  samples <- tstrsplit(uni_id_column, sep)[[1]]
  cells   <- tstrsplit(uni_id_column, sep)[[2]]
  fmt_ids <- paste(cells,samples,sep = "-")
  
  return(fmt_ids)
}

bcr_df$unique_cell_id <- formatUniqueIDcolum(bcr_df$unique_cell_id)


# Generate Hamming distances
bcr_results_igh <- bcr_df[bcr_df$locus %in% "IGH",]
dist_nearest    <- distToNearest(bcr_results_igh)

# Automated Threshold - will take a long time with the full cohort
# - Threshold setting occurs by default when running the immcantation pipeline
threshold_output <- shazam::findThreshold(dist_nearest$dist_nearest,
                                          method = "gmm", model = "gamma-norm",
                                          cutoff = "user",spc = 0.99)
threshold <- threshold_output@threshold

p_ham_auto <- plot(threshold_output, binwidth = 0.02, silent = TRUE) + theme(axis.title = element_text(size = 18))


ham_plot_file_name <- paste("full_cohort_distn_hamming_CDR3_auto_thresh_",threshold,".pdf",sep = "")
pdf(file = file.path(plot_dir,ham_plot_file_name),height = 5,width = 7)
plot(p_ham_auto)
dev.off()

