#!/usr/bin/env Rscript

# START ----
# Workflow for creating analysis objects
#  - Prior to running this code set working directory to 'sc-cart-analysis/code'

library(Seurat)
library(rjson)
library(tidyverse)
library(stringr)
library(knitr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(kableExtra)
library(devtools)
library(parallel)
library(data.table)
library(harmony)
library(DoubletFinder)
library(impactSingleCellToolkit)

# INPUT and OUTPUT Directories
param_file_fh = "../input/input_full_cohort_analysis_cart.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads

# OUTPUT Results Directories
plot_dir     = params$plot_directory
objects_dir  = params$objects_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}

# RData Objects to save
fh_raw_seurat_obj <- file.path(objects_dir,"seurat_obj.combined.RData")
fh_raw_csp_mtx <- file.path(objects_dir,"csp_count_matrix.combined.RData")

# Increase size of ENV
options(future.globals.maxSize= 891289600)


# 1. Separate samples by batch/assay ----

# create batch:[gex/csf/multi] column
metadata$run_group <- paste(metadata$run,metadata$type,sep = ":")
metadata_gex       <- metadata[metadata$type %in% c("counts","counts_gex"),]
metadata_csp       <- metadata[metadata$type %in% c("counts","counts_csp"),]

run_metadata_gex_list <- split(metadata_gex,metadata_gex$run_group)
run_metadata_csp_list <- split(metadata_csp,metadata_csp$run_group)


# 2. Create Seurat Objects for each batch ----

# makeRunInputMtx
# - output.type = c("gex","csp")
makeRunInputMtx <- function(runID ,run_metadata_list,THREADS,
                            output.type = "gex",min.genes.gex=700,min.genes.csp=9) {
  run          <- tstrsplit(runID,":")[[1]]
  type         <- tstrsplit(runID,":")[[2]]
  multi.status <- FALSE
  if(type == "counts"){
    multi.status <- TRUE
  }
  
  run_metadata <- run_metadata_list[[runID]]
  run_metadata <- run_metadata[run_metadata$type %in% type,]
  
  # Get run directory and sample vector
  sample_dir_vec <- run_metadata$results_directory_path
  dataset_loc    <- tstrsplit(sample_dir_vec,"multi_counts/")[[1]] %>% unique()
  dataset_loc    <- paste(dataset_loc,"multi_counts/",sep = "")
  
  # Separate samples based on whether they are mult, gex only or csp only
  samples.vec   <- run_metadata$sample
  
  if(THREADS > length(samples.vec)){
    THREADS <- length(samples.vec)
  }
  min_genes <- 0
  if(output.type == "gex"){min_genes <- min.genes.gex}
  if(output.type == "csp"){min_genes <- min.genes.csp}
  gex.matrix <- generateCombinedMatrix(dataset_loc, samples.vec,THREADS = THREADS,multi.results = multi.status,
                                       assay = output.type,min.genes.per.cell = min_genes,max.genes.per.cell = NULL)
  
  return(gex.matrix)
}

gex_run_list <- names(run_metadata_gex_list) %>% as.list()
csp_run_list <- names(run_metadata_csp_list) %>% as.list()

# Initialize variables
dataset_loc        <- ""
samples.vec        <- c()
multi.results      <- NULL
assay              <- NULL
min.genes.per.cell <- NULL
max.genes.per.cell <- NULL

# CSP input matrix
csp_mtx_list <- lapply(csp_run_list, makeRunInputMtx,
                       run_metadata_list=run_metadata_csp_list,
                       THREADS=THREADS,output.type = "csp",min.genes.csp=1)
merged_csp_mtx <- do.call("cbind",csp_mtx_list)

# GEX input matrix
gex_mtx_list <- lapply(gex_run_list, makeRunInputMtx,
                       run_metadata_list=run_metadata_gex_list,
                       THREADS=THREADS,output.type = "gex",min.genes.gex=400)
merged_gex_mtx <- do.call("cbind",gex_mtx_list)



# 3. Create Seurat Objects for GEX and CSP ----
seurat.obj <- createSeuratObjectFromMatrix(
  sc.data      = merged_gex_mtx,
  project.name = "GEX_CART",
  npca         = 20, min.genes = 400,
  normalize = F,dim.reduction = F
)

# Add samples column
cells      <- row.names(seurat.obj@meta.data)
sample_col <- cells
sample_col <- tstrsplit(sample_col,"-")[[2]]
seurat.obj@meta.data$sample <- sample_col

# Omit samples with low number of cells
LOW_CELLS_THRESH          <- 100
sample_cell_counts        <- as.data.frame(table(seurat.obj$sample))
names(sample_cell_counts) <- c("sample","Freq")
sample_cell_counts        <- sample_cell_counts[sample_cell_counts$Freq > LOW_CELLS_THRESH,]

cells_to_keep  <- names(seurat.obj$sample[seurat.obj$sample %in% sample_cell_counts$sample])
seurat.obj     <- seurat.obj[,colnames(seurat.obj) %in% cells_to_keep]

# 4. Identify Doublets ----
set.seed(1234)

seurat.obj <- findDoublets(seurat.obj,sample.col = "sample",threads = THREADS)

# 5. Add in CSP and overlap both assays ----
cells_csp <- colnames(merged_csp_mtx)
cells_gex <- colnames(seurat.obj)
overlap   <- cells_csp[cells_csp %in% cells_gex]

# Only keep overlapping cells
seurat.obj <- seurat.obj[,colnames(seurat.obj) %in% overlap]
counts_csp <- merged_csp_mtx[,colnames(merged_csp_mtx) %in% overlap]

seurat.obj[["CSP"]] <- CreateAssayObject(counts = counts_csp)


# SAVE ----
save(seurat.obj,file = fh_raw_seurat_obj)
save(merged_csp_mtx,file = fh_raw_csp_mtx)
