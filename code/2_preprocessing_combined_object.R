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

# Input and Output RData Objects
fh_raw_seurat_obj       <- file.path(objects_dir,"seurat_obj.combined.RData")
fh_processed_seurat_obj <- file.path(objects_dir,"seurat_obj.processed.rds")

# Increase size of ENV
options(future.globals.maxSize= 40000*1024^2)
set.seed(112358)



# 1. Load and start normalization/clustering ----

# Load combined GEX/CSP object
load(fh_raw_seurat_obj)

# Omit TCR V genes from object
all_genes  <- row.names(seurat.obj)
tcr_vgenes <- all_genes[grep("TRAV|TRBV|TRGV|TRDV",all_genes)]
seurat.obj <- seurat.obj[! row.names(seurat.obj) %in% tcr_vgenes,]

# Normalization per sample
if(PERSAMPLE_SCT){
  seurat.obj.list <- SplitObject(seurat.obj, split.by="sample")
  seurat.obj.list <- lapply(seurat.obj.list, function(x) {
    SCTransform(x, return.only.var.genes = FALSE, verbose = FALSE)
  })
  
  features <- SelectIntegrationFeatures(object.list = seurat.obj.list, nfeatures = 2000)
  
  # Set variable features for each object to avoid SCT model conflicts
  for (i in 1:length(seurat.obj.list)) {
    VariableFeatures(seurat.obj.list[[i]]) <- features
  }
  
  # Step 5: Prepare for integration
  seurat.obj.list <- PrepSCTIntegration(object.list = seurat.obj.list, anchor.features = features)
  
  # Find integration anchors
  reference_samples <- NULL
  #if(anchor_col %in% names(metadata)){
  #  reference_samples        <- 1:length(names(seurat.obj.list))
  #  names(reference_samples) <- names(seurat.obj.list)
  #  
  #  status_col   <- metadata[,anchor_col,drop=T]
  #  meta_anchors <- metadata[ status_col %in% "yes",]
  #  meta_anchors <- meta_anchors[meta_anchors$type %in% c("counts","counts_gex"),]
  #  
  #  reference_samples        <- reference_samples[names(reference_samples) %in% meta_anchors$sample] %>% as.numeric()
  #}
  
  anchors <- FindIntegrationAnchors(object.list = seurat.obj.list,
                                    reference = reference_samples,
                                    normalization.method = "SCT", 
                                    anchor.features = features)
  
  # Integrate data
  # - integration layer will only have the variable to genes, but other layers will have all genes
  # - To include all genes, features.to.integrate = all_genes # all_genes is a vector will all gene names
  seurat.obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  #DefaultAssay(seurat.obj) <- "integrated"
}


# Re-cluster
norm_param   <- TRUE
if(PERSAMPLE_SCT){norm_param   <- FALSE}

seurat.obj   <- scaleAndClusterSeuratObject(seurat.obj,normalize = norm_param,dims = 1:30,npca = 10,tsne = T)
pc_elbowplot <- plotOptimalPCsforSeuratObject(seurat.obj)

### Temporarily save this object
save(seurat.obj, file = file.path(objects_dir,"seurat_obj.TEMP.RData"))

# 2. Plot pre-batch correction and add in batch column ----

# add batch column
metadata_gex                <- metadata[metadata$type %in% c("counts_gex","counts"),]
sample_to_batch_list        <- metadata_gex$run
names(sample_to_batch_list) <- metadata_gex$sample
sample_to_batch_list        <- as.list(sample_to_batch_list)

batch_col   <- seurat.obj$sample
uni_samples <- unique(batch_col)
for (s in uni_samples) {
  batch <- sample_to_batch_list[[s]]
  batch_col[batch_col %in% s] <- batch
}
seurat.obj$batch <- batch_col

# Plot 
sample_umap <- UMAPPlot(object = seurat.obj, group.by="sample")

png(file.path(plot_dir,"sample_umap.png"),height = 500,width = 600)
print(sample_umap)
dev.off()


# 3. Azimuth automated celltype annotations ----
###InstallData("pbmcref")
###DefaultAssay(seurat.obj). # Check that assay is either RNA or SCT, not CSP

# Run Azimuth
results_azimuth  <- RunAzimuth(seurat.obj, reference = "pbmcref")

save(results_azimuth, file = file.path(objects_dir,"results_azimuth.TEMP.RData"))

# Add annotations back into Seurat object
celltypes_azimuth <- results_azimuth@meta.data
annot_cols        <- names(celltypes_azimuth)
annot_cols        <- annot_cols[annot_cols %in% c("predicted.celltype.l1","predicted.celltype.l2","predicted.celltype.l3")]
celltypes_azimuth <- celltypes_azimuth[,annot_cols]

seurat.obj@meta.data <- cbind(seurat.obj@meta.data,celltypes_azimuth)


# 4. UMAP celltype annotations ----

# Plot Azimuth annotations

### Level 1 annotations (main)
l1_cell_annot_azimuth <- UMAPPlot(object = seurat.obj, group.by="predicted.celltype.l1")

png(file.path(plot_dir,"celltype_annot_azimuth_l1_umap.png"),height = 500,width = 600)
print(l1_cell_annot_azimuth)
dev.off()

### Level 2 annotations (fine)
l2_cell_annot_azimuth <- UMAPPlot(object = seurat.obj, group.by="predicted.celltype.l2")

png(file.path(plot_dir,"celltype_annot_azimuth_l2_umap.png"),height = 500,width = 800)
print(l2_cell_annot_azimuth)
dev.off()


# 6. Process CSP part of object ----
# the CSP layer was removed during per sample integration

#DefaultAssay(seurat.obj) <- "CSP" # it should be 'SCT' after processing the GEX data

### Normalize
#seurat.obj <- NormalizeData(seurat.obj, normalization.method = "CLR", margin = 2, assay = "CSP")
### Run PCA
#VariableFeatures(seurat.obj) = rownames(seurat.obj@assays[["CSP"]])
#seurat.obj = seurat.obj %>% 
#  ScaleData(verbose=F) %>%
#  RunPCA(reduction.name="apca",approx=F, verbose=F) 
### Pick PCs
#total_variance <- sum(matrixStats::rowVars(
#  as.matrix(seurat.obj@assays[["CSP"]]@scale.data)))
#eigValues = (seurat.obj@reductions$apca@stdev)^2  
#varExplained = eigValues / total_variance
#csp_pc_plot <- plot(varExplained)
### Run UMAP
#seurat.obj <- RunUMAP(seurat.obj, 
#              reduction = 'apca', 
#              dims = 1:12, 
#              assay = 'CSP', 
#              reduction.name = 'csp.umap', 
#              reduction.key = 'cspUMAP_',
#              verbose = F)

# SAVE ----
DefaultAssay(seurat.obj) <- "SCT"
saveRDS(seurat.obj,file = fh_processed_seurat_obj)
