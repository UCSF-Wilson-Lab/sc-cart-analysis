#!/usr/bin/env Rscript

# START ----
# Workflow for processing combined analysis object
#  - Prior to running this code set working directory to 'sc-cart-analysis/code'

library(Seurat)
library(SeuratData)
library(rjson)
library(tidyverse)
library(stringr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(Azimuth)
library(kableExtra)
library(devtools)
library(parallel)
library(data.table)
library(harmony)
library(DoubletFinder)
library(impactSingleCellToolkit)

# INPUT and OUTPUT Directories
# 2_preprocessing_combined_object_gex.R  [Param File]

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("ERROR: At least one argument must be supplied (JSON parameter file).json", call.=FALSE)
} 
param_file_fh = args[1]
#param_file_fh = "../input/input_full_cohort_analysis_cart.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads
PERSAMPLE_SCT = TRUE

# OUTPUT Results Directories
plot_dir     = params$plot_directory
objects_dir  = params$objects_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}

# Input and Output RData/RDS Objects
input_csp_mtx_fh        <- file.path(objects_dir,"csp_count_matrix.combined.RData")
fh_processed_seurat_obj <- file.path(objects_dir,"seurat_obj.processed.rds")

# 1. Load Objects
load(input_csp_mtx_fh)
seurat.obj <- readRDS(fh_processed_seurat_obj)


# 2. Add in CSP and overlap both assays ----
cells_csp <- colnames(merged_csp_mtx)
cells_gex <- colnames(seurat.obj)
overlap   <- cells_csp[cells_csp %in% cells_gex]

# Only keep overlapping cells
seurat.obj <- seurat.obj[,colnames(seurat.obj) %in% overlap]
counts_csp <- merged_csp_mtx[,colnames(merged_csp_mtx) %in% overlap]

seurat.obj[["CSP"]] <- CreateAssayObject(counts = counts_csp)

# 3. Process CSP part of object ----
DefaultAssay(seurat.obj) <- "CSP" # it should be 'SCT' after processing the GEX data

### Normalize
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "CLR", margin = 2, assay = "CSP")
### Run PCA
VariableFeatures(seurat.obj) = rownames(seurat.obj@assays[["CSP"]])
seurat.obj = seurat.obj %>% 
  ScaleData(verbose=F) %>%
  RunPCA(reduction.name="apca",approx=F, verbose=F) 
### Pick PCs
total_variance <- sum(matrixStats::rowVars(
  as.matrix(seurat.obj@assays[["CSP"]]@scale.data)))
eigValues = (seurat.obj@reductions$apca@stdev)^2  
varExplained = eigValues / total_variance
csp_pc_plot <- plot(varExplained)
### Run UMAP
seurat.obj <- RunUMAP(seurat.obj, 
              reduction = 'apca', 
              dims = 1:12, 
              assay = 'CSP', 
              reduction.name = 'csp.umap', 
              reduction.key = 'cspUMAP_',
              verbose = F)

# SAVE ----
DefaultAssay(seurat.obj) <- "SCT"
saveRDS(seurat.obj,file = fh_processed_seurat_obj)
