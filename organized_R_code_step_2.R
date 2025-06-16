# NOTES --------------------------------------------------------------------------
# Script: Single-Cell RNA-seq Analysis in normal and unloaded (BTX) tendon enthesis groups @ p7, p14, and p28
# the named rds file can be found in /rds file/normal_vs_unloaded_entheses_@_p7_p14_p28_annotated.rds
# Organization Date: 06-16-2025
# Description: This script performs loading, QC, doublet removal (DoubletFinder),
#              ambient RNA correction (SoupX), batch correction (scVI),
#              clustering, and cell type annotation for scRNA-seq data from
#              10X Genomics and BGI platforms.
# ALL CELL CLUSTERS--------------------------------------------------------------------------

# --- 0.1 Environment Setup ---
# Set working directory (IMPORTANT: Modify this to your actual base directory)
# project_base_dir <- ''
# setwd(project_base_dir)

rm(list = ls())
gc() 

# --- 0.2 Load Libraries ---
# Core R packages
library(Matrix)
library(magrittr) # For pipe operator %>%
library(dplyr)
library(stringr)
library(tibble)
library(tidyverse) # Includes dplyr, tibble, stringr, ggplot2 etc.
library(data.table)

# Seurat and related scRNA-seq packages
library(Seurat)
library(SeuratDisk) # For h5ad conversion
library(scCustomize) # Custom Seurat plotting (optional, if used for specific plots)

# Specific analysis packages
library(DoubletFinder)
library(SoupX)

# Reticulate for Python integration (scVI)
library(reticulate)
library(sceasy)    # For Seurat to AnnData conversion
library(anndata)   # Used by sceasy

# Plotting
library(ggsci)     # For color palettes (e.g., pal_simpsons)
library(RColorBrewer) # For additional color palettes

# --- 0.3 Configuration & Global Parameters ---
# Directories (adjust BASE_DIR if not using setwd())
BASE_DIR <- if (exists("project_base_dir")) project_base_dir else "~/compare/"
DATA_DIR <- file.path(BASE_DIR, "data")
RESULTS_DIR <- file.path(BASE_DIR, "results_organized") # Store outputs here
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
  message(paste("Created results directory:", RESULTS_DIR))
}

# QC Parameters
MIN_FEATURE_RNA_QUANTILE <- 0.05
MAX_PERCENT_MT <- 10

# Gene patterns/lists for filtering (from expression matrix)
MT_GENE_PATTERN <- "^mt-"
RIBO_GENE_PATTERN <- "^Rp[sl]" # For mouse, human is "^RP[SL]"
HSP_GENE_PATTERN <- "^Hsp"    # Case-sensitive, adjust if needed
DISSOCIATION_GENES <- c(
  'Atf3','Zfp36','Fos','Dusp1','Egr1','Gadd45g','Jun','Junb','Jund',
  'Actg1','Ier2','Sat1','Ubc','Nr4a1','Egr2','Fosb','Ppp1r15a',
  'Mt1','Mt2','Errfi1','Ier3','Nfkbia','Zfp36l1','Cebpb','Cebpd','Btg2'
) %>% unique() # Ensure no duplicates in the list

# DoubletFinder Parameters
DOUBLET_RATE <- 0.076 # Expected doublet rate
DF_ASSAY <- "SCT"     # Assay for DoubletFinder preprocessing
DF_PCS <- 1:30        # PCs for DoubletFinder (original used 1:50 after SCT, adjust as needed)
DF_RESOLUTION <- 0.4  # Clustering resolution for homotypic proportion model

# scVI Parameters (IMPORTANT: Ensure this Python environment and path are correct)
SCVI_PYTHON_PATH <- Sys.getenv("RETICULATE_PYTHON", unset="/media/whc/02C401E8C401DF33/zt/.local/share/r-miniconda/envs/scvi/bin/python")
SCVI_CONDA_ENV <- "scvi"
SCVI_BATCH_KEY <- "sample_id" # Metadata column for scVI batch correction (will be created)

# Clustering Parameters
FINAL_CLUSTER_RESOLUTION <- 0.3 # Resolution for final clustering post-scVI

# File Names for RDS objects and other outputs
BEFORE_QC_RDS <- file.path(RESULTS_DIR, "seurat_merged_before_QC.rds")
QC_LIST_RDS <- file.path(RESULTS_DIR, "seurat_list_after_QC.rds")
DOUBLETFINDER_LIST_RDS <- file.path(RESULTS_DIR, "seurat_list_after_doubletfinder.rds")
SOUPX_DETAILS_RDS <- file.path(RESULTS_DIR, "soupX_channel_object.rds")
SOUPX_CORRECTED_RDS <- file.path(RESULTS_DIR, "seurat_merged_after_soupX.rds")
SCVI_INTEGRATED_RDS <- file.path(RESULTS_DIR, "seurat_merged_after_scVI.rds")
FINAL_ANNOTATED_RDS <- file.path(RESULTS_DIR, "seurat_merged_final_annotated.rds")
ALL_MARKERS_CSV <- file.path(RESULTS_DIR, paste0("all_markers_res", FINAL_CLUSTER_RESOLUTION, ".csv"))

# Output gene list files (cumulative unique genes removed across samples)
MT_GENES_FILE <- file.path(RESULTS_DIR, "cumulative_mitochondrial_genes_removed.csv")
RIBO_GENES_FILE <- file.path(RESULTS_DIR, "cumulative_ribosomal_genes_removed.csv")
HSP_GENES_FILE <- file.path(RESULTS_DIR, "cumulative_hsp_genes_removed.csv")
DISS_GENES_FILE <- file.path(RESULTS_DIR, "cumulative_dissociation_genes_removed.csv")

# Clean up previous gene list files if they exist to avoid appending to old data
file.remove(MT_GENES_FILE, RIBO_GENES_FILE, HSP_GENES_FILE, DISS_GENES_FILE, showWarnings = FALSE)

# --- 0.4 Helper Functions ---

#' Load 10X Genomics or BGI scRNA-seq data.
#' Assumes folder names are used as project identifiers.
load_scrna_data <- function(data_dir, folder_pattern, data_type_prefix, gene_col_for_bgi = 1) {
  folders <- list.files(data_dir, pattern = folder_pattern, full.names = FALSE)
  if (length(folders) == 0) {
    warning(paste("No folders found matching pattern:", folder_pattern, "in directory:", data_dir))
    return(list())
  }
  
  sce_list <- lapply(folders, function(folder_name) {
    message(paste("Processing folder:", file.path(data_dir, folder_name)))
    counts_path <- file.path(data_dir, folder_name)
    
    current_counts <- Read10X(counts_path, 
                              gene.column = if (grepl("bgi", folder_pattern, ignore.case = TRUE)) gene_col_for_bgi else 2)
    
    seurat_obj <- CreateSeuratObject(counts = current_counts, project = folder_name) # folder_name becomes orig.ident
    return(seurat_obj)
  })
  names(sce_list) <- folders # Name list elements for easier access
  return(sce_list)
}

#' Filter specified genes from a Seurat object's RNA assay and record them.
filter_and_record_genes <- function(seurat_obj, gene_pattern = NULL, gene_list = NULL, output_file, gene_set_name) {
  rna_assay <- seurat_obj@assays$RNA
  if (is.null(rna_assay)) {
    message(sprintf("RNA assay not found for sample %s. Skipping %s gene filtering.", 
                    unique(seurat_obj$orig.ident), gene_set_name))
    return(seurat_obj)
  }
  all_genes <- rownames(rna_assay)
  
  genes_to_remove_indices <- c()
  if (!is.null(gene_pattern)) {
    genes_to_remove_indices <- c(genes_to_remove_indices, grep(gene_pattern, all_genes, ignore.case = TRUE))
  }
  if (!is.null(gene_list)) {
    genes_to_remove_indices <- c(genes_to_remove_indices, which(all_genes %in% gene_list))
  }
  genes_to_remove_indices <- unique(genes_to_remove_indices)
  
  if (length(genes_to_remove_indices) > 0) {
    genes_identified <- all_genes[genes_to_remove_indices]
    
    # Append unique genes found in this sample to the global list (if not already there)
    if (length(genes_identified) > 0) {
      if (!file.exists(output_file)) {
        write.table(data.frame(Gene = character()), output_file, sep = ",", col.names = TRUE, row.names = FALSE, quote=FALSE)
      }
      existing_genes <- tryCatch(read.csv(output_file, stringsAsFactors = FALSE)$Gene, error = function(e) character(0))
      new_genes_to_write <- setdiff(genes_identified, existing_genes)
      if(length(new_genes_to_write) > 0) {
        write.table(data.frame(Gene = new_genes_to_write), output_file, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE, quote=FALSE)
      }
    }
    
    features_to_keep <- setdiff(all_genes, genes_identified)
    # Subset the RNA assay directly if only counts and data slots are relevant
    if (length(features_to_keep) < length(all_genes)) {
      seurat_obj@assays$RNA <- subset(rna_assay, features = features_to_keep)
      message(sprintf("Removed %d %s genes from RNA assay of sample %s. %d features remaining in RNA assay.",
                      length(genes_identified), gene_set_name, unique(seurat_obj$orig.ident), length(features_to_keep)))
    }
  } else {
    message(sprintf("No %s genes to remove found in sample %s.", gene_set_name, unique(seurat_obj$orig.ident)))
  }
  return(seurat_obj)
}

# --- 1. Data Loading and Initial Merging ---
message("--- Part 1: Data Loading and Initial Merging ---")

# 1.1 Load 10X Genomics data
message("Loading 10X Genomics data...")
sce_list_10x <- load_scrna_data(DATA_DIR, folder_pattern = "10x_filtered_feature_bc_matrix", data_type_prefix = "10X")

# 1.2 Load BGI data
message("Loading BGI data...")
sce_list_bgi <- load_scrna_data(DATA_DIR, folder_pattern = "bgi_filtered_feature_bc_matrix", data_type_prefix = "BGI", gene_col_for_bgi = 1)

# 1.3 Combine lists and merge
all_sce_list <- c(sce_list_10x, sce_list_bgi)
if (length(all_sce_list) == 0) {
  stop("No data loaded. Please check DATA_DIR and folder_patterns.")
}

message("Merging all datasets...")
# Add cell IDs prefixes to ensure unique cell names after merging
# Project names (orig.ident) are already folder names from load_scrna_data
if (length(all_sce_list) == 1) {
  seurat_merged_raw <- all_sce_list[[1]]
} else {
  seurat_merged_raw <- merge(all_sce_list[[1]], 
                             y = all_sce_list[2:length(all_sce_list)],
                             add.cell.ids = names(all_sce_list)) # Use folder names as prefixes
}
seurat_merged_raw@project.name <- "Compare_Raw_Merged"

message(paste("Saving merged raw data to:", BEFORE_QC_RDS))
saveRDS(seurat_merged_raw, BEFORE_QC_RDS)

# 1.4 Split merged object by original identity for individual QC
message("Splitting merged object for individual QC...")
sce_list_for_qc <- SplitObject(seurat_merged_raw, split.by = 'orig.ident')

rm(sce_list_10x, sce_list_bgi, all_sce_list, seurat_merged_raw); gc()


# --- 2. Quality Control (per sample) ---
message("--- Part 2: Quality Control (per sample) ---")
# Optional: Load if starting from here
# sce_list_for_qc <- SplitObject(readRDS(BEFORE_QC_RDS), split.by = 'orig.ident')

sce_list_filtered <- lapply(X = sce_list_for_qc, FUN = function(x) {
  current_ident <- as.character(x$orig.ident[1]) # orig.ident is the folder name
  message(paste("Performing QC for sample:", current_ident))
  
  # 2.1 Add metadata based on orig.ident (folder name)
  # Assumes folder names like "btx1w_10x_filtered_feature_bc_matrix" or "nom2w_bgi_..."
  x$sample_id <- substr(current_ident, 1, 5) # e.g., "btx1w"
  x$group_id  <- substr(current_ident, 1, 3) # e.g., "btx"
  x$time_id   <- substr(current_ident, 4, 5) # e.g., "1w"
  x$tech_id   <- substr(current_ident, 7, 9) # e.g., "10x" or "bgi"
  
  # 2.2 Calculate mitochondrial and ribosomal content
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = MT_GENE_PATTERN, assay = "RNA")
  x[["percent.ribo"]] <- PercentageFeatureSet(x, pattern = RIBO_GENE_PATTERN, assay = "RNA")
  
  # 2.3 Filter cells based on nFeature_RNA and percent.mt
  nfeature_rna_low_threshold <- quantile(x$nFeature_RNA, probs = MIN_FEATURE_RNA_QUANTILE, na.rm = TRUE)
  cells_before_filter <- ncol(x)
  x <- subset(x, subset = nFeature_RNA > nfeature_rna_low_threshold & percent.mt < MAX_PERCENT_MT)
  message(sprintf("Filtered cells for %s: %d to %d (nFeature_RNA > %.0f, percent.mt < %d)",
                  current_ident, cells_before_filter, ncol(x), nfeature_rna_low_threshold, MAX_PERCENT_MT))
  
  if (ncol(x) == 0) {
    message(paste("Warning: All cells filtered out for sample", current_ident, "- skipping gene filtering for this sample."))
    return(NULL) # Return NULL if no cells remain
  }
  
  # 2.4 Filter specific gene sets from the RNA assay's expression matrix
  # This removes the genes themselves, not just cells expressing them.
  x <- filter_and_record_genes(x, gene_pattern = MT_GENE_PATTERN, output_file = MT_GENES_FILE, gene_set_name = "Mitochondrial")
  x <- filter_and_record_genes(x, gene_pattern = RIBO_GENE_PATTERN, output_file = RIBO_GENES_FILE, gene_set_name = "Ribosomal")
  x <- filter_and_record_genes(x, gene_pattern = HSP_GENE_PATTERN, output_file = HSP_GENES_FILE, gene_set_name = "Heat Shock Protein")
  x <- filter_and_record_genes(x, gene_list = DISSOCIATION_GENES, output_file = DISS_GENES_FILE, gene_set_name = "Dissociation-induced")
  
  return(x)
})

# Remove NULL elements from the list (samples with no cells after QC)
sce_list_filtered <- sce_list_filtered[!sapply(sce_list_filtered, is.null)]
if(length(sce_list_filtered) == 0) {
  stop("All samples were filtered out during QC. Check QC parameters and input data.")
}

message(paste("Saving QC-filtered Seurat list to:", QC_LIST_RDS))
saveRDS(sce_list_filtered, QC_LIST_RDS)
rm(sce_list_for_qc); gc()


# --- 3. Doublet Detection and Removal (per sample) ---
message("--- Part 3: Doublet Detection and Removal (per sample) ---")
# Optional: Load if starting from here
# sce_list_filtered <- readRDS(QC_LIST_RDS)

sce_list_after_doubletfinder <- lapply(X = sce_list_filtered, FUN = function(x) {
  current_ident <- unique(x$orig.ident)
  message(paste("Running DoubletFinder for sample:", current_ident))
  
  # Preprocessing for DoubletFinder (SCTransform, PCA, UMAP, Clustering)
  # variable.features.n = 3000-4000 is common. Original had 4000.
  x <- SCTransform(x, assay="RNA", new.assay.name = "SCT", vst.flavor = "v2", variable.features.n = 4000, verbose = FALSE)
  
  x <- RunPCA(x, assay = DF_ASSAY, npcs = max(DF_PCS), verbose = FALSE)
  x <- RunUMAP(x, assay = DF_ASSAY, dims = DF_PCS, reduction.name = "umap.DF", reduction.key = "UMAPDF_", verbose = FALSE) # Separate UMAP for DF
  x <- FindNeighbors(x, assay = DF_ASSAY, dims = DF_PCS, verbose = FALSE)
  x <- FindClusters(x, resolution = DF_RESOLUTION, verbose = FALSE) # Used for homotypic proportion
  
  # Optimize pK value
  sweep.res.list <- paramSweep_v3(x, PCs = DF_PCS, sct = TRUE) # sct=TRUE as data is SCTransformed
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  
  optimal_pK <- 0.09 # Default pK
  if (nrow(bcmvn) > 0 && any(!is.na(bcmvn$BCmetric))) {
    optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  } else {
    message(paste("pK optimization failed for", current_ident, "- using default pK =", optimal_pK))
  }
  
  # Estimate expected number of doublets
  annotations <- x$seurat_clusters # Using seurat_clusters as proxy for cell types
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(DOUBLET_RATE * ncol(x)) 
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  x <- doubletFinder_v3(x, PCs = DF_PCS, pN = 0.25, pK = optimal_pK, 
                        nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  
  # Identify the classification column (DoubletFinder names it based on parameters)
  df_classification_col <- colnames(x@meta.data)[grep(paste0("^DF.classifications_0.25_", optimal_pK), colnames(x@meta.data))]
  if (length(df_classification_col) == 0) {
    # Fallback to general pattern if specific PK not found (e.g. if PK was rounded differently in colname)
    df_classification_col <- colnames(x@meta.data)[grep("^DF.classifications", colnames(x@meta.data))]
    df_classification_col <- df_classification_col[length(df_classification_col)] # Take the last one
  }
  if (length(df_classification_col) != 1) {
    stop(paste("Could not uniquely identify DoubletFinder classification column for", current_ident))
  }
  
  # Subset to keep singlets
  x$doublet_classification <- x@meta.data[[df_classification_col]]
  singlets <- subset(x, subset = doublet_classification == "Singlet")
  message(sprintf("Identified doublets for %s. Kept %d singlets out of %d cells.",
                  current_ident, ncol(singlets), ncol(x)))
  
  # Return only the RNA assay with original (pre-SCT) counts for singlets
  # This is important for SoupX which needs raw counts.
  # SCT data is not needed beyond this point for this sample.
  rna_assay_singlets <- singlets@assays$RNA # This assay still holds original counts
  meta_data_singlets <- singlets@meta.data
  
  # Create a new Seurat object with only the RNA assay of singlets
  # This ensures clean data for SoupX.
  # Make sure the counts are the raw counts, not normalized data.
  # `SCTransform` does not overwrite `RNA@counts` by default.
  singlet_obj_rna_only <- CreateSeuratObject(counts = rna_assay_singlets@counts,
                                             meta.data = meta_data_singlets,
                                             project = unique(singlets$orig.ident))
  return(singlet_obj_rna_only)
})

# Remove NULLs again if any sample processing failed or resulted in no singlets
sce_list_after_doubletfinder <- sce_list_after_doubletfinder[!sapply(sce_list_after_doubletfinder, is.null)]
if(length(sce_list_after_doubletfinder) == 0) {
  stop("All samples were filtered out after DoubletFinder. Check parameters.")
}

message(paste("Saving Seurat list (RNA assay of singlets) after DoubletFinder to:", DOUBLETFINDER_LIST_RDS))
saveRDS(sce_list_after_doubletfinder, DOUBLETFINDER_LIST_RDS)
rm(sce_list_filtered, sweep.res.list, sweep.stats, bcmvn); gc()


# --- 4. Ambient RNA Correction with SoupX ---
message("--- Part 4: Ambient RNA Correction with SoupX ---")
# sce_list_after_doubletfinder <- readRDS(DOUBLETFINDER_LIST_RDS)
# seurat_merged_raw <- readRDS(BEFORE_QC_RDS) # Contains all droplets

# 4.1 Merge the QCed, doublet-filtered list for SoupX input (filtered cell counts = toc)
message("Merging doublet-filtered samples for SoupX 'toc' matrix...")
if (length(sce_list_after_doubletfinder) == 1) {
  seurat_for_toc <- sce_list_after_doubletfinder[[1]]
} else {
  # Need to get unique names for add.cell.ids if list elements are not named
  # sce_list_after_doubletfinder should be named by orig.ident from previous step
  # If not, generate unique IDs
  sample_names_for_merge <- sapply(sce_list_after_doubletfinder, function(s) unique(s$orig.ident))
  if (any(duplicated(sample_names_for_merge))) { # Should not happen if orig.ident is unique
    sample_names_for_merge <- make.unique(sample_names_for_merge)
  }
  seurat_for_toc <- merge(sce_list_after_doubletfinder[[1]], 
                          y = sce_list_after_doubletfinder[2:length(sce_list_after_doubletfinder)],
                          add.cell.ids = sample_names_for_merge) # Use unique sample names for cell IDs
}
toc <- seurat_for_toc@assays$RNA@counts # Table Of Counts (filtered cells)
message(sprintf("Dimensions of 'toc' (post-QC, post-doublet cell counts): %d genes x %d cells", nrow(toc), ncol(toc)))

# 4.2 Get the full count matrix (before QC) for SoupX input (all droplet counts = tod)
message("Loading raw merged data for SoupX 'tod' matrix...")
seurat_merged_raw <- readRDS(BEFORE_QC_RDS)
DefaultAssay(seurat_merged_raw) <- "RNA"
tod <- seurat_merged_raw@assays$RNA@counts   # Table Of Droplets (all barcodes)
message(sprintf("Dimensions of 'tod' (raw droplet counts): %d genes x %d cells", nrow(tod), ncol(tod)))

# Ensure genes in tod and toc match. SoupX expects common genes.
# Genes were already filtered in QC step from 'toc' (e.g. mt, ribo).
# So, 'tod' must be subsetted to these genes.
common_genes <- intersect(rownames(tod), rownames(toc))
if (length(common_genes) == 0) stop("No common genes between tod and toc for SoupX.")
tod <- tod[common_genes, ]
toc <- toc[common_genes, ] # toc should already have these genes, but good to ensure
message(sprintf("Dimensions of 'tod' after common gene selection: %d x %d", nrow(tod), ncol(tod)))
message(sprintf("Dimensions of 'toc' after common gene selection: %d x %d", nrow(toc), ncol(toc)))

cells_passed_qc_and_doublet <- colnames(toc) # These are prefixed names if merge was used.

common_cells <- intersect(colnames(tod), colnames(toc))
if (length(common_cells) != ncol(toc)) {
  message(sprintf("Warning: Mismatch in cell numbers. %d cells in toc, %d common with tod. Using common cells.", ncol(toc), length(common_cells)))
  toc <- toc[, common_cells]
  seurat_for_toc <- subset(seurat_for_toc, cells = common_cells) # subset seurat object too
  if (ncol(toc) == 0) stop("No common cells between tod and toc after name reconciliation.")
}

# 4.3 Estimate contamination using SoupX
message("Performing quick clustering on 'toc' for SoupX cell groups...")
# Create a temporary Seurat object from 'toc' for quick clustering.
# `seurat_for_toc` already has the counts for `toc` and its metadata.
quick_cluster_obj <- seurat_for_toc # Use the Seurat object that `toc` came from
DefaultAssay(quick_cluster_obj) <- "RNA"
quick_cluster_obj <- NormalizeData(quick_cluster_obj, verbose = FALSE)
quick_cluster_obj <- FindVariableFeatures(quick_cluster_obj, nfeatures = 2000, verbose = FALSE)
quick_cluster_obj <- ScaleData(quick_cluster_obj, features = VariableFeatures(quick_cluster_obj), verbose = FALSE)
quick_cluster_obj <- RunPCA(quick_cluster_obj, features = VariableFeatures(quick_cluster_obj), npcs = 20, verbose = FALSE)
quick_cluster_obj <- FindNeighbors(quick_cluster_obj, dims = 1:15, verbose = FALSE)
quick_cluster_obj <- FindClusters(quick_cluster_obj, resolution = 0.3, verbose = FALSE)

soup_clusters <- quick_cluster_obj@meta.data$seurat_clusters
names(soup_clusters) <- rownames(quick_cluster_obj@meta.data) # names are cell barcodes

message("Initializing SoupChannel and estimating contamination...")
sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE) # Profile helps autoEstCont
sc <- setClusters(sc, soup_clusters) # Provide cell clusters
sc <- autoEstCont(sc, verbose = TRUE) # Automatically estimate contamination
estimated_rho <- sc$metaData$rho
message(paste("SoupX estimated contamination fraction (rho):", round(estimated_rho, 3)))

message(paste("Saving SoupX channel object to:", SOUPX_DETAILS_RDS))

message("Adjusting counts with SoupX...")
# `toc` (and thus `seurat_for_toc`) cell names should match `sc$toc` cell names
out_counts <- adjustCounts(sc, roundToInt = TRUE) # roundToInt for downstream analysis
message(sprintf("Dimensions of SoupX adjusted counts: %d genes x %d cells", nrow(out_counts), ncol(out_counts)))

# 4.4 Create Seurat object with corrected counts
# Metadata should come from `seurat_for_toc` (which matches `toc` cells)
seurat_corrected <- CreateSeuratObject(counts = out_counts, meta.data = seurat_for_toc@meta.data)
seurat_corrected@project.name <- "Compare_SoupX_Corrected"

# Factor re-leveling for metadata (sample_id, time_id, group_id were added in QC)
# Ensure these levels match your actual data
expected_sample_levels <- c('btx1w','nom1w','btx2w','nom2w','btx4w','nom4w')
expected_time_levels <- c('1w','2w','4w')
expected_group_levels <- c('nom','btx')

seurat_corrected$sample_id <- factor(x = seurat_corrected$sample_id, levels = intersect(expected_sample_levels, unique(seurat_corrected$sample_id)))
seurat_corrected$time_id <- factor(x = seurat_corrected$time_id, levels = intersect(expected_time_levels, unique(seurat_corrected$time_id)))
seurat_corrected$group_id <- factor(x = seurat_corrected$group_id, levels = intersect(expected_group_levels, unique(seurat_corrected$group_id)))

message(paste("Saving SoupX-corrected Seurat object to:", SOUPX_CORRECTED_RDS))
saveRDS(seurat_corrected, SOUPX_CORRECTED_RDS)
rm(seurat_merged_raw, seurat_for_toc, tod, toc, quick_cluster_obj, sc, out_counts, sce_list_after_doubletfinder); gc()


# --- 5. Batch Correction and Dimensionality Reduction with scVI ---
message("--- Part 5: Batch Correction with scVI ---")
# Optional: Load if starting here
# seurat_corrected <- readRDS(SOUPX_CORRECTED_RDS)

# 5.1 Setup Python environment for scVI
message(paste("Using Python for scVI from conda env:", SCVI_CONDA_ENV, "or path:", SCVI_PYTHON_PATH))
# Check if RETICULATE_PYTHON is set, otherwise use SCVI_PYTHON_PATH
# `use_python` or `use_condaenv` should be called *before* importing Python modules.
if (Sys.getenv("RETICULATE_PYTHON") == "") {
  use_python(SCVI_PYTHON_PATH, required = TRUE)
} else {
  use_condaenv(SCVI_CONDA_ENV, required = TRUE) # This is often more robust
}
py_config() # Display Python configuration

scvi <- import("scvi", convert = FALSE)
sc <- import("scanpy", convert = FALSE)

# 5.2 Convert Seurat object to AnnData
message("Converting Seurat object to AnnData format...")
# scVI works best on raw/adjusted counts. `seurat_corrected` has SoupX-adjusted counts in RNA assay.
adata <- sceasy::convertFormat(seurat_corrected, from = "seurat", to = "anndata",
                               assay = "RNA", main_layer = "counts", # Use "counts" layer
                               drop_single_values = FALSE)
print(adata)

# 5.3 Setup AnnData for scVI model
# `SCVI_BATCH_KEY` (e.g., "sample_id") must be a column in `adata.obs` (from `seurat_corrected@meta.data`)
message(paste("Setting up AnnData for scVI model with batch key:", SCVI_BATCH_KEY))
if (!SCVI_BATCH_KEY %in% colnames(seurat_corrected@meta.data)) {
  stop(paste("Batch key '", SCVI_BATCH_KEY, "' not found in Seurat object metadata.", sep=""))
}
scvi$model$SCVI$setup_anndata(adata, batch_key = SCVI_BATCH_KEY)

# 5.4 Create and train scVI model
message("Creating and training scVI model (this may take time)...")
model <- scvi$model$SCVI(adata)
model$train() # Default training parameters
message("scVI model training complete.")

# 5.5 Extract latent representation
message("Extracting latent representation from scVI model...")
latent_representation <- model$get_latent_representation()
latent_matrix <- as.matrix(py_to_r(latent_representation))
rownames(latent_matrix) <- colnames(seurat_corrected) # Assign cell names

# 5.6 Add scVI embeddings to Seurat object
seurat_integrated <- seurat_corrected
seurat_integrated[["scvi"]] <- CreateDimReducObject(embeddings = latent_matrix,
                                                    key = "scvi_", # Embeddings will be scvi_1, scvi_2, ...
                                                    assay = DefaultAssay(seurat_integrated)) # Associated with RNA assay

message(paste("Saving Seurat object with scVI embeddings to:", SCVI_INTEGRATED_RDS))
saveRDS(seurat_integrated, SCVI_INTEGRATED_RDS)
rm(seurat_corrected, adata, model, latent_representation, latent_matrix); gc()


# --- 6. Clustering, Visualization, and Annotation (Post-scVI) ---
message("--- Part 6: Clustering, Visualization, and Annotation ---")

# 6.1 SCTransform for downstream analysis (e.g. marker finding, visualization)
# Run SCTransform on the SoupX-corrected counts ("RNA" assay).
message("Running SCTransform on SoupX-corrected counts for marker analysis and visualization...")
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- SCTransform(seurat_integrated, vst.flavor = "v2", verbose = FALSE,
                                 assay = "RNA", new.assay.name = "SCT_SOUPX") # Store in new assay
# DefaultAssay(seurat_integrated) <- "SCT_SOUPX" # Set as default for subsequent steps like marker finding

# 6.2 Clustering using scVI embeddings
message("Performing clustering using scVI embeddings...")
scvi_dims <- 1:ncol(Embeddings(seurat_integrated, reduction = "scvi"))

message("Running UMAP on scVI embeddings...")
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "scvi", dims = scvi_dims,
                             reduction.name = "umap.scvi", reduction.key = "UMAPscvi_")

message("Finding neighbors and clusters on scVI embeddings...")
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "scvi", dims = scvi_dims)
seurat_integrated <- FindClusters(seurat_integrated, resolution = FINAL_CLUSTER_RESOLUTION, verbose = FALSE)
# Cluster IDs are in seurat_integrated$seurat_clusters by default

# Visualization of clusters
umap_plot_clusters <- DimPlot(seurat_integrated, reduction = "umap.scvi", group.by = "seurat_clusters",
                              label = TRUE, label.size = 5, pt.size = 0.5, repel = TRUE) +
  ggtitle(paste("UMAP of scVI Integrated Data (Resolution", FINAL_CLUSTER_RESOLUTION, ")"))
print(umap_plot_clusters)
ggsave(file.path(RESULTS_DIR, paste0("umap_scvi_clusters_res", FINAL_CLUSTER_RESOLUTION, ".png")), plot = umap_plot_clusters, width = 10, height = 8)

# 6.3 Marker Gene Identification
message("Identifying cluster markers using SCT_SOUPX assay...")
DefaultAssay(seurat_integrated) <- "SCT_SOUPX" # Ensure this assay is used for markers
Idents(seurat_integrated) <- "seurat_clusters" # Set idents to the new clusters
seurat_integrated <- PrepSCTFindMarkers(seurat_integrated, assay = "SCT_SOUPX", verbose = FALSE)

all_markers <- FindAllMarkers(seurat_integrated,
                              assay = "SCT_SOUPX", # Specify assay
                              only.pos = TRUE,
                              max.cells.per.ident = 200, # As in original
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              verbose = FALSE)

# Save markers
if (nrow(all_markers) > 0) {
  write.csv(all_markers, ALL_MARKERS_CSV, row.names = FALSE)
  message(paste("Saved all cluster markers to:", ALL_MARKERS_CSV))
  # Print top markers for a few clusters
  message("Top markers per cluster (example):")
  print(all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))
} else {
  message("No markers found with current parameters.")
}

# Example FeaturePlots and VlnPlots (using genes from original script)
genes_to_plot <- c('Col10a1','Tnn','Scx','Tnmd','Ly6a','Sox9','Lepr','Runx2','Sfrp5')
genes_present_in_data <- genes_to_plot[genes_to_plot %in% rownames(seurat_integrated[["SCT_SOUPX"]])]

if (length(genes_present_in_data) > 0) {
  message(paste("Plotting features for:", paste(genes_present_in_data, collapse=", ")))
  fp <- FeaturePlot(seurat_integrated, features = genes_present_in_data, reduction = "umap.scvi", order = TRUE, pt.size=0.2)
  print(fp)
  ggsave(file.path(RESULTS_DIR, "featureplot_selected_genes.png"), plot = fp, width = 12, height = max(8, 2*ceiling(length(genes_present_in_data)/3)))
  
  if ("Tnn" %in% genes_present_in_data && "time_id" %in% colnames(seurat_integrated@meta.data)) {
    vp <- VlnPlot(seurat_integrated, features = "Tnn", split.by = "time_id", group.by="seurat_clusters",
                  assay = "SCT_SOUPX", pt.size = 0) +
      ggtitle("Tnn Expression by Time Point and Cluster")
    print(vp)
    ggsave(file.path(RESULTS_DIR, "vlnplot_Tnn_by_time_cluster.png"), plot = vp, width = 10, height = 6)
  }
} else {
  message("None of the specified genes for FeaturePlot are present in the SCT_SOUPX assay.")
}

# 6.4 Cell Type Annotation (example, adapt based on your markers)
message("Performing cell type annotation (example renaming)...")
# This renaming map is highly specific to your dataset and clustering outcome.
# It will likely need adjustment if `FINAL_CLUSTER_RESOLUTION` or data changes.
# Original cluster numbers: '0'-'17'
rename_map <- c(
  `0`='Articular chondrocyte (intermediate zone)', 
  `1`='Growth plate chondrocyte_1', # _1 to distinguish from other GP chond.
  `2`='Tendon and related cell',
  `3`='B cell_1', 
  `4`='B cell_2', 
  `5`='Articular chondrocyte (superficial zone)', 
  `6`='Growth plate chondrocyte_2',
  `7`='B cell_3',
  `8`='B cell_4',
  `9`='Bone marrow cell',
  `10`='Erythroid-like cell_1', 
  `11`='Macrophage',
  `12`='B cell_5', 
  `13`='Smooth muscle cell', 
  `14`='NK and T cell',
  `15`='Endothelial cell',
  `16`='Erythroid-like cell_2',
  `17`='Neural cell'
)

# Apply renaming if clusters exist
current_cluster_ids <- levels(seurat_integrated$seurat_clusters)
valid_rename_map <- rename_map[names(rename_map) %in% current_cluster_ids]
missing_from_map <- setdiff(current_cluster_ids, names(rename_map))
if(length(missing_from_map) > 0) {
  message(paste("Warning: Clusters", paste(missing_from_map, collapse=", "), "are not in the rename_map and will keep their numeric ID."))
  # Add them to map to keep original ID
  names(missing_from_map) <- missing_from_map
  valid_rename_map <- c(valid_rename_map, missing_from_map)
}

seurat_integrated <- RenameIdents(seurat_integrated, valid_rename_map)
seurat_integrated$cell_type_level1 <- Idents(seurat_integrated) # Store in metadata

# Visualization with new cell type names
num_cell_types <- length(levels(seurat_integrated$cell_type_level1))
cell_type_colors <- if (num_cell_types <= 16) {
  pal_simpsons("springfield")(num_cell_types) # From ggsci
} else {
  colorRampPalette(RColorBrewer::brewer.pal(n = max(3, min(12, num_cell_types)), name = "Paired"))(num_cell_types)
}
names(cell_type_colors) <- levels(seurat_integrated$cell_type_level1)


umap_plot_celltypes <- DimPlot(seurat_integrated, reduction = "umap.scvi", group.by = "cell_type_level1",
                               label = TRUE, repel = TRUE, label.size = 3.5, pt.size = 0.5,
                               cols = cell_type_colors) +
  ggtitle("UMAP with Annotated Cell Types") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size=3))) # Improve legend
print(umap_plot_celltypes)
ggsave(file.path(RESULTS_DIR, "umap_scvi_annotated_celltypes.png"), plot = umap_plot_celltypes, width = 12, height = 10)


# --- 7. Final Save and Cleanup ---
message("--- Part 7: Final Save and Cleanup ---")

message(paste("Saving final annotated Seurat object to:", FINAL_ANNOTATED_RDS))
saveRDS(seurat_integrated, FINAL_ANNOTATED_RDS)

# Final garbage collection
gc()

sessionInfo()
message(paste0("Script execution finished successfully at ", Sys.time(), "! Results are in: ", RESULTS_DIR))



# TENDON AND ENTHESIS DOWNSTREAM ANALYSIS--------------------------------------------------------------------------
# Script: Focused Analysis of Tendon and Enthesis Cells
# Organization Date: 06-16-2025
# Description: This script takes a subset of cells (Tendon and Enthesis related)
#              from a larger Seurat object, re-integrates them using scVI,
#              performs detailed clustering, visualization (DimPlots, FeaturePlots, SCP plots),
#              trajectory analysis (Slingshot via SCP), GO enrichment,
#              and other specialized analyses like CytoTRACE2, miloR, ROGUE, and ClusterGVis.
# Prerequisite: A Seurat object named /rds file/normal_vs_unloaded_entheses_@_p7_p14_p28_annotated.rds
#               with a metadata column 'name1' containing cell type annotations.


compare_all <- readRDS('./normal_vs_unloaded_entheses_@_p7_p14_p28_annotated.rds')
# --- 8.1 Subset Cells for Focused Analysis ---
message("--- Part 8.1: Subsetting 'Tendon and enthesis cell' ---")

if (exists("compare_all") && "name1" %in% colnames(compare_all@meta.data)) {
  Idents(compare_all) <- 'name1'
  if (SUBSET_CELL_TYPE_NAME %in% levels(Idents(compare_all))) {
    related_initial <- subset(compare_all, idents = SUBSET_CELL_TYPE_NAME)

    if ("RNA" %in% Assays(related_initial)) {
      related <- CreateSeuratObject(counts = related_initial@assays$RNA@counts,
                                    meta.data = related_initial@meta.data)
      related@project.name <- "Tendon_Enthesis_Subset"
      message(sprintf("Subset '%s' created with %d cells and %d features.",
                      SUBSET_CELL_TYPE_NAME, ncol(related), nrow(related)))
      # The line `related$ <- NULL` seems like a typo or an incomplete thought. Removed.
      
      saveRDS(related, RELATED_INITIAL_SUBSET_RDS)
    } else {
      stop("RNA assay not found in the subset. Cannot proceed with scVI.")
    }
  } else {
    stop(paste("Cell type '", SUBSET_CELL_TYPE_NAME, "' not found in 'name1' idents of 'compare_all'.", sep=""))
  }
  rm(related_initial); gc()
} else {
  stop("'compare_all' object or 'name1' metadata column not found. Cannot perform subsetting.")
}

# --- 8.2 Re-integrate Subset with scVI ---
message("--- Part 8.2: Re-integrating subset with scVI ---")
DefaultAssay(related) <- 'RNA' # Ensure RNA assay (with counts) is default

# Setup Python environment for scVI
message(paste("Using Python for scVI from conda env:", SCVI_CONDA_ENV_SUBSET, "or path:", SCVI_PYTHON_PATH_SUBSET))
# Sys.unsetenv("RETICULATE_PYTHON") # This line might interfere if RETICULATE_PYTHON was set globally for a reason. Use with caution.
use_python(SCVI_PYTHON_PATH_SUBSET, required = TRUE) # Or use_condaenv(SCVI_CONDA_ENV_SUBSET, required = TRUE)
py_config()

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

# Convert Seurat object to AnnData
adata <- sceasy::convertFormat(related, from="seurat", to="anndata",
                               assay="RNA", main_layer="counts", drop_single_values=FALSE)
print(adata)

# Setup AnnData for scVI model (no batch key specified in original code for subset)
scvi$model$SCVI$setup_anndata(adata) # If batch correction is needed for subset, add batch_key here

# Create and train scVI model
model <- scvi$model$SCVI(adata)
message(paste("Training scVI model for subset with max_epochs:", SCVI_MAX_EPOCHS_SUBSET))
model$train(max_epochs = as.integer(SCVI_MAX_EPOCHS_SUBSET))

# Get latent representation
latent <- model$get_latent_representation()
latent_matrix <- as.matrix(py_to_r(latent))
rownames(latent_matrix) <- colnames(related)

# Add scVI embeddings to Seurat object
related[["scvi_subset"]] <- CreateDimReducObject(embeddings = latent_matrix, key = "scvisub_", assay = DefaultAssay(related))
message("scVI integration for subset complete.")
rm(adata, model, latent, latent_matrix); gc()

# --- 8.3 Clustering and Visualization of scVI-integrated Subset ---
message("--- Part 8.3: Clustering and Visualization of scVI-integrated Subset ---")
scvi_subset_dims <- 1:ncol(Embeddings(related, reduction = "scvi_subset"))

related <- FindNeighbors(related, dims = scvi_subset_dims, reduction = "scvi_subset", graph.name = "nn.scvi_subset")
related <- FindClusters(related, resolution = SUBSET_CLUSTER_RESOLUTION, graph.name = "nn.scvi_subset") # Store clusters in related$seurat_clusters

related <- RunUMAP(related, dims = scvi_subset_dims, reduction = "scvi_subset",
                   n.components = 2, reduction.name = "umap.scvi_subset", reduction.key = "UMAPsub_")

# Plot UMAP
subset_umap_plot <- DimPlot(related, reduction = "umap.scvi_subset", group.by = "seurat_clusters",
                            label = TRUE, label.size = 5, pt.size = 1, repel = TRUE) +
  ggtitle(paste("UMAP of Subset (scVI, Res:", SUBSET_CLUSTER_RESOLUTION, ")"))
print(subset_umap_plot)
ggsave(file.path(RESULTS_DIR_SUBSET, paste0("umap_subset_scvi_res", SUBSET_CLUSTER_RESOLUTION, ".png")),
       plot = subset_umap_plot, width = 8, height = 7)


DefaultAssay(related) <- "RNA"
related <- NormalizeData(related, verbose = FALSE)
related <- FindVariableFeatures(related, verbose = FALSE)
related <- ScaleData(related, features = VariableFeatures(related), verbose = FALSE)


Idents(related) <- "seurat_clusters"
subset_markers_file <- file.path(RESULTS_DIR_SUBSET, paste0("markers_subset_res", SUBSET_CLUSTER_RESOLUTION, ".csv"))
if (!file.exists(subset_markers_file)) { # Avoid re-running if file exists
  # Assuming RNA assay was LogNormalized and Scaled above:
  subset_markers <- FindAllMarkers(related, assay = "RNA", only.pos = TRUE, max.cells.per.ident = 200, logfc.threshold = 0.25, min.pct = 0.1)
  write.csv(subset_markers, subset_markers_file, row.names = FALSE)
  message(paste("Subset cluster markers saved to:", subset_markers_file))
} else {
  subset_markers <- read.csv(subset_markers_file)
  message(paste("Loaded subset cluster markers from:", subset_markers_file))
}


# FeaturePlot examples
genes_to_plot_subset <- c('Scx','Sox9','Col2a1','Tnn','Tnmd','Ly6a','Clec3a','Wif1','Mkx')
genes_present_subset <- genes_to_plot_subset[genes_to_plot_subset %in% rownames(related)]
if (length(genes_present_subset) > 0) {
  fp_subset <- FeaturePlot(related, features = genes_present_subset, reduction = "umap.scvi_subset", order = TRUE, pt.size=0.5)
  print(fp_subset)
  ggsave(file.path(RESULTS_DIR_SUBSET, "featureplot_subset_selected_genes.png"), plot = fp_subset, width = 12, height = max(8, 2*ceiling(length(genes_present_subset)/3)))
}

saveRDS(related, RELATED_SCVI_INTEGRATED_RDS)
# related <- readRDS(RELATED_SCVI_INTEGRATED_RDS) # Optional load


# Re-factor to ensure correct levels after subsetting.
if ("sample_id" %in% colnames(related@meta.data)) {
  related$sample_id <- factor(related$sample_id, levels = intersect(c('btx1w','nom1w','btx2w','nom2w','btx4w','nom4w'), unique(related$sample_id)))
}
if ("group_id" %in% colnames(related@meta.data)) {
  related$group_id <- factor(related$group_id, levels = intersect(c('btx','nom'), unique(related$group_id)))
}

saveRDS(related, RELATED_FINAL_SUBSET_RDS)
# related <- readRDS(RELATED_FINAL_SUBSET_RDS) # Optional load


# --- 8.5 Label Prediction and Visualization (scanvi/external labels) ---
message("--- Part 8.5: Label Prediction/Loading and Visualization ---")
# This section assumes an AnnData file with predictions (e.g., from scanVI)
# or that 'related' already has a 'predictions_scanvi' column from a loaded RDS.

if (file.exists(RELATED_PREDICTED_ADATA_PATH)) {
  message(paste("Loading predicted labels from AnnData file:", RELATED_PREDICTED_ADATA_PATH))
  # Ensure reticulate is configured
  # use_python(SCVI_PYTHON_PATH_SUBSET, required = TRUE) # If not already done
  
  # Convert AnnData to Seurat
  # The output file name in convertFormat is just for the temporary RDS, not the object name.
  related_predicted_temp_obj <- sceasy::convertFormat(RELATED_PREDICTED_ADATA_PATH, from="anndata", to="seurat",
                                                      outFile = RELATED_PREDICTED_RDS_TEMP)
  
  # Check if 'predictions_scanvi' exists in the loaded object's metadata
  if ("predictions_scanvi" %in% colnames(related_predicted_temp_obj@meta.data)) {
    # Align cells: 'related' (current subset) vs 'related_predicted_temp_obj'
    common_cells_pred <- intersect(colnames(related), colnames(related_predicted_temp_obj))
    if (length(common_cells_pred) > 0) {
      related <- subset(related, cells = common_cells_pred) # Subset current object
      predictions_df <- related_predicted_temp_obj@meta.data[common_cells_pred, "predictions_scanvi", drop = FALSE]
      related <- AddMetaData(related, metadata = predictions_df)
      related$name2 <- related$predictions_scanvi # Rename column
      message("Added 'name2' predictions to 'related' object.")
    } else {
      warning("No common cells between current 'related' object and the loaded predictions. Skipping label transfer.")
    }
  } else {
    warning(paste("'predictions_scanvi' column not found in loaded AnnData from", RELATED_PREDICTED_ADATA_PATH))
  }
  rm(related_predicted_temp_obj); file.remove(RELATED_PREDICTED_RDS_TEMP)
} else {
  # Fallback or alternative: check if 'related' was loaded with 'name2' (e.g. from 'compare_related_CCA_RNA_scalded_named.rds')
  if (!"name2" %in% colnames(related@meta.data)) {
    message(paste("AnnData file for predictions not found at:", RELATED_PREDICTED_ADATA_PATH, "and 'name2' not in metadata. Skipping label prediction part or manual annotation needed."))
  } else {
    message("'name2' column already exists in 'related' object metadata.")
  }
}

# If 'name2' column is now present:
if ("name2" %in% colnames(related@meta.data)) {
  related$name2 <- factor(x=related$name2,
                          levels = intersect(c('Enthesis chondrocytes','Tenocytes','related progenitors','Tenoblasts','Mesenchymal progenitors'), unique(related$name2)))
  
  # Normalization and Scaling (if not done already or if RNA assay needs reprocessing for 'name2' markers)
  DefaultAssay(related) <- 'RNA'
  related <- NormalizeData(related, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  related <- FindVariableFeatures(related, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  related <- ScaleData(related, features = rownames(related), verbose = TRUE) # Simpler scaling
  
  saveRDS(related, RELATED_WITH_PREDICTIONS_RDS)
  # related <- readRDS(RELATED_WITH_PREDICTIONS_RDS)
  
  # Define colors for 'name2'
  name2_levels <- levels(related$name2)
  name2_cols <- ggsci::pal_locuszoom("default")(length(name2_levels))
  names(name2_cols) <- name2_levels
  
  # DimPlot using scCustomize
  if (length(name2_levels) > 0) {
    pdf(file = file.path(RESULTS_DIR_SUBSET, 'umap_subset_name2_scCustom.pdf'), width = 7.5, height = 4.5)
    print(DimPlot_scCustom(seurat_object = related, colors_use = name2_cols, reduction = "umap.scvi_subset",
                           group.by = "name2", label = TRUE, label.size = 4, pt.size = 1.5))
    dev.off()
    
    # CellDimPlot from scplotter (part of SCP)
    pdf(file = file.path(RESULTS_DIR_SUBSET, 'umap_subset_name2_SCP.pdf'), width = 8, height = 8)
    print(scplotter::CellDimPlot(related, group_by = "name2", reduction = "umap.scvi_subset",
                                 pt_size = 1.2, label = TRUE, label_size = 3, label_fg = "black", label_bg = "white",
                                 palcolor = name2_cols, theme_light = FALSE)) 
    dev.off()
    
    # CellStatPlot for cell numbers
    pdf(file = file.path(RESULTS_DIR_SUBSET, 'cellstat_name2_by_sample_SCP.pdf'), width = 7, height = 5)
    print(SCP::CellStatPlot(related, stat.by = "name2", group.by = "sample_id", palcolor = name2_cols, label = TRUE,
                            base_size = 10, aspect.ratio = 1.5))
    dev.off()
  }
  
  # FeatureDimPlot (SCP)
  genes_for_scp_fp <- c("Sox9", "Scx", "Gli1", "Clec3a",'Ly6a','Prrx1')
  genes_present_scp_fp <- genes_for_scp_fp[genes_for_scp_fp %in% rownames(related)]
  if (length(genes_present_scp_fp) > 0) {
    pdf(file = file.path(RESULTS_DIR_SUBSET, 'featureplot_subset_name2_genes_SCP.pdf'), width = 9, height = 6)
    print(SCP::FeatureDimPlot(srt = related, features = genes_present_scp_fp, ncol = 3, reduction = "umap.scvi_subset", theme_light = FALSE))
    dev.off()
  }
  
  # FindAllMarkers for 'name2'
  Idents(related) <- 'name2'
  DefaultAssay(related) <- 'RNA' # Use RNA for interpretability if scaled appropriately
  name2_markers_file <- file.path(RESULTS_DIR_SUBSET, "markers_subset_name2.csv")
  if(!file.exists(name2_markers_file)){
    related_name2_markers <- FindAllMarkers(related, assay = "RNA", only.pos = TRUE,
                                            min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(related_name2_markers, name2_markers_file, row.names=FALSE)
  } else {
    related_name2_markers <- read.csv(name2_markers_file)
  }
} else {
  message("Skipping visualizations and marker finding for 'name2' as the column is not available.")
}

# --- 8.6 Trajectory Analysis with Slingshot (via SCP) ---
message("--- Part 8.6: Trajectory Analysis with Slingshot ---")
if ("name2" %in% colnames(related@meta.data) && "umap.scvi_subset" %in% Reductions(related)) {

  start_node_trajectory <- 'Mesenchymal progenitors'
  end_nodes_trajectory <- c('Enthesis chondrocytes', 'Tenocytes')
  
  # Check if start/end nodes are present in 'name2' levels
  if (start_node_trajectory %in% levels(related$name2) && all(end_nodes_trajectory %in% levels(related$name2))) {
    message("Running Slingshot for trajectory inference...")
    related <- SCP::RunSlingshot(srt = related, group.by = "name2",
                                 reduction = "umap.scvi_subset", # Use the UMAP from subset scVI
                                 start_node = start_node_trajectory,
                                 end_node = end_nodes_trajectory,
                                 show_plot = FALSE) # Set to TRUE to see default plot during run
    
    # Visualize lineages
    # Identify lineage names generated by Slingshot (e.g., Lineage1, Lineage2)
    lineage_names <- names(related@tools)[grep("^Slingshot$", names(related@tools))][[1]] # Get the tool name
    lineage_names_actual <- names(slot(related@tools[[lineage_names]], "lineages"))
    
    if (length(lineage_names_actual) > 0) {
      # Define colors for lineages dynamically or use fixed ones if known
      lineage_colors <- RColorBrewer::brewer.pal(n = max(3, length(lineage_names_actual)), name = "Set1")[1:length(lineage_names_actual)]
      
      pdf(file = file.path(RESULTS_DIR_SUBSET, 'trajectory_name2_lineages_SCP.pdf'), width = 8, height = 7)
      print(scplotter::CellDimPlot(related, group_by = "name2", reduction = "umap.scvi_subset",
                                   lineages = lineage_names_actual,
                                   lineages_palcolor = lineage_colors, # Auto or specify
                                   palcolor = name2_cols, pt_size =1.1,
                                   theme_light = FALSE))
      dev.off()
      
      # Dynamic Gene Features (example for first two lineages)
      DefaultAssay(related) <- 'RNA' # Ensure RNA assay for expression
      lineages_for_dynamic_genes <- lineage_names_actual[1:min(2, length(lineage_names_actual))]
      if (length(lineages_for_dynamic_genes) > 0) {
        related <- SCP::RunDynamicFeatures(srt = related, lineages = lineages_for_dynamic_genes, n_candidates = 200) # Original had 1000, then 200
        
        # Plot dynamic genes (example)
        dynamic_genes_to_plot <- c('Tnmd',"Wif1",'Scrg1','Tnn') # Example from original
        dynamic_genes_present <- dynamic_genes_to_plot[dynamic_genes_to_plot %in% rownames(related)]
        
        if (length(dynamic_genes_present) > 0) {
          pdf(file = file.path(RESULTS_DIR_SUBSET, 'trajectory_name2_dynamic_genes_SCP.pdf'), width = 6, height = 3 * length(dynamic_genes_present))
          print(SCP::DynamicPlot(srt = related, lineages = lineages_for_dynamic_genes, group.by = "name2",
                                 line_palcolor = lineage_colors[1:length(lineages_for_dynamic_genes)],
                                 point_palcolor = name2_cols,
                                 features = dynamic_genes_present, ncol = 1,
                                 compare_lineages = TRUE, compare_features = FALSE))
          dev.off()
        }
      }
    } else {
      message("Slingshot did not identify any lineages.")
    }
    saveRDS(related, RELATED_TRAJECTORY_RDS)
    # related <- readRDS(RELATED_TRAJECTORY_RDS)
  } else {
    message("Start or end nodes for trajectory not found in 'name2' levels. Skipping Slingshot.")
  }
} else {
  message("Skipping trajectory analysis: 'name2' column or 'umap.scvi_subset' reduction not found.")
}


# --- 8.7 Specialized Analyses (CytoTRACE2, miloR, ROGUE, etc.) ---
message("--- Part 8.7: Specialized Analyses ---")

# 8.7.1 CytoTRACE2 (Potency analysis)
if ("RNA" %in% Assays(related)) {
  message("Running CytoTRACE2...")
  tryCatch({
    # if (!requireNamespace("CytoTRACE2", quietly = TRUE)) devtools::install_github("digitalCytometry/CytoTRACE2")
    library(CytoTRACE2)
    related_ct <- cytotrace2(related, assay_ epigenomics = "RNA", # Original code used 'is_seurat = TRUE'
                             slot_type = "counts", species = 'mouse', seed =  1234)
    related$CytoTRACE2_Score <- related_ct$CytoTRACE2_Score
    related$CytoTRACE2_Potency <- related_ct$CytoTRACE2_Potency # Potency category
    
    pdf(file.path(RESULTS_DIR_SUBSET, "cytotrace2_plots.pdf"), width=7, height=5)
    print(SCP::FeatureStatPlot(related, stat.by = "CytoTRACE2_Score", group.by = "sample_id", plot_type = "violin") + ggtitle("CytoTRACE2 Score by Sample"))
    print(SCP::CellStatPlot(related, stat.by = "CytoTRACE2_Potency", group.by = "sample_id", label = TRUE) + ggtitle("CytoTRACE2 Potency by Sample"))
    dev.off()
    saveRDS(related, RELATED_CYTOTRACE_RDS)
    # related <- readRDS(RELATED_CYTOTRACE_RDS)
  }, error = function(e) {
    message(paste("CytoTRACE2 failed:", e$message))
  })
}

# 8.7.2 miloR (Differential Abundance)
if ("scvi_subset" %in% Reductions(related) && "sample_id" %in% colnames(related@meta.data) && "group_id" %in% colnames(related@meta.data) && "time_id" %in% colnames(related@meta.data)) {
  message("Running miloR for differential abundance...")
  tryCatch({
    library(miloR)
    library(SingleCellExperiment)
    library(scater) # for plotReducedDim
    
    related_sce <- as.SingleCellExperiment(related, assay = "RNA") # Use RNA counts
    colData(related_sce)$PCA <- Embeddings(related, "scvi_subset") # Use scVI latent space as PCA for Milo
    
    milo_obj <- Milo(related_sce)
    milo_obj <- buildGraph(milo_obj, k = 11, d = ncol(Embeddings(related, "scvi_subset")), reduced.dim = "PCA") # d is number of scVI dims
    milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 11, d = ncol(Embeddings(related, "scvi_subset")), refined = TRUE, reduced.dim = "PCA")
    
    pdf(file.path(RESULTS_DIR_SUBSET, "miloR_nhood_size_hist.pdf"), width=6, height=4)
    print(plotNhoodSizeHist(milo_obj))
    dev.off()
    
    milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), sample="sample_id")
    
    milo_design <- data.frame(colData(milo_obj))[,c("sample_id", "group_id",'time_id')] %>% distinct()
    rownames(milo_design) <- milo_design$sample_id
    # Ensure factors are set correctly for the model
    milo_design$group_id <- factor(milo_design$group_id)
    milo_design$time_id <- factor(milo_design$time_id)
    
    # This needs to be adapted if your sample order or number is different.
    # For now, using as is:
    if (nrow(milo_design) >=6) { # Check if enough samples for this specific slicing
      milo_design_ordered <- milo_design[match(c('btx1w','nom1w','btx2w','nom2w','btx4w','nom4w'), milo_design$sample_id),]
      milo_design_ordered <- milo_design_ordered[!is.na(milo_design_ordered$sample_id),] # Keep only existing samples
    } else {
      milo_design_ordered <- milo_design
    }
    
    da_results <- testNhoods(milo_obj, design = ~ group_id + time_id, # Simplified model, adjust as needed
                             design.df = milo_design_ordered)
    
    pdf(file.path(RESULTS_DIR_SUBSET, "miloR_pvalue_hist.pdf"), width=6, height=4)
    print(ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) + ggtitle("MiloR P-value Distribution"))
    dev.off()
    
    milo_obj <- buildNhoodGraph(milo_obj)
    colData(milo_obj)$UMAP <- Embeddings(related, "umap.scvi_subset") # Add UMAP for plotting
    
    # Plot Milo results on UMAP
    pdf(file.path(RESULTS_DIR_SUBSET, 'miloR_DA_on_UMAP.pdf'), width = 8, height = 7)
    print(plotNhoodGraphDA(milo_obj, da_results, layout="UMAP", alpha=0.1) + 
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
            ggtitle("Differential Abundance (MiloR)"))
    dev.off()
    
    # Beeswarm plot (if 'name2' is available for annotation)
    if ("name2" %in% colnames(related@meta.data)) {
      da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "name2")
      pdf(file.path(RESULTS_DIR_SUBSET, 'miloR_DA_beeswarm.pdf'), width = 7, height = 6)
      print(ggplot(da_results, aes(logFC, name2_fraction, color = logFC)) + 
              ggbeeswarm::geom_quasirandom(size =2) +
              scale_color_gradient2(low = "blue", mid = 'grey', high = "red", midpoint = 0)+
              theme_bw() + labs(title="DA by Dominant Cell Type", x="logFC (group_id effect)", y = "Dominant Cell Type in Nhood")
      )
      dev.off()
    }
    saveRDS(milo_obj, RELATED_MILO_OBJECT_RDS)
    saveRDS(da_results, RELATED_MILO_RESULTS_RDS)
  }, error = function(e) {
    message(paste("miloR failed:", e$message))
  })
}

# 8.7.3 ROGUE (Cell Purity/Heterogeneity)
if ("RNA" %in% Assays(related) && "name2" %in% colnames(related@meta.data) && "sample_id" %in% colnames(related@meta.data)) {
  message("Running ROGUE for cell type purity...")
  tryCatch({
    library(ROGUE)
    rogue_input_matrix <- related@assays$RNA@counts 
    
    # Ensure labels and samples are character vectors
    rogue_labels <- as.character(related$name2)
    rogue_samples <- as.character(related$sample_id)
    
    # Check for NAs and remove corresponding cells
    valid_indices <- !is.na(rogue_labels) & !is.na(rogue_samples)
    rogue_input_matrix_valid <- rogue_input_matrix[, valid_indices]
    rogue_labels_valid <- rogue_labels[valid_indices]
    rogue_samples_valid <- rogue_samples[valid_indices]
    
    if (ncol(rogue_input_matrix_valid) > 0 && length(unique(rogue_labels_valid)) > 1 && length(unique(rogue_samples_valid)) > 0) {
      rogue.res <- rogue(rogue_input_matrix_valid, labels = rogue_labels_valid,
                         samples = rogue_samples_valid, platform = "UMI", span = 0.6) 
      
      pdf(file.path(RESULTS_DIR_SUBSET, 'rogue_purity_boxplot.pdf'), height = 5, width = max(6, length(unique(rogue_samples_valid))*1.5))
      print(rogue.boxplot(rogue.res) + theme_bw() + ylab("ROGUE Score") +
              theme(axis.text.x = element_text(angle=45, hjust=1)))
      dev.off()
    } else {
      message("Skipping ROGUE due to insufficient valid data (NAs, too few unique labels or samples).")
    }
  }, error = function(e) {
    message(paste("ROGUE failed:", e$message))
  })
}

# 8.7.4 ClusterGVis (Gene Set Enrichment and Visualization)
if (exists("related_name2_markers") && nrow(related_name2_markers) > 0 && "name2" %in% colnames(related@meta.data)) {
  message("Running ClusterGVis for enrichment visualization...")
  tryCatch({
    library(ClusterGVis)
    library(org.Mm.eg.db) # For mouse
    
    # Prepare data for ClusterGVis (using top N markers per cluster for 'name2')
    top_markers_for_gvis <- related_name2_markers %>%
      dplyr::filter(p_val_adj < 0.05) %>% # Filter significant markers
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 15, wt = avg_log2FC) # Take top 15 for example
    
    if (nrow(top_markers_for_gvis) > 0) {
      # Use SCT assay if available and scaled, otherwise RNA (scaled)
      assay_for_gvis <- if ("SCT_subset" %in% Assays(related)) "SCT_subset" else "RNA"
      
      st.data <- prepareDataFromscRNA(object = related,
                                      group.by = "name2",
                                      assays = assay_for_gvis, # Make sure this assay has scaled data
                                      diffData = top_markers_for_gvis,
                                      keep.uniqGene = FALSE, # Keep all provided markers
                                      showAverage = TRUE) # Use average expression
      
      # Perform enrichment
      enrich.res <- enrichCluster(object = st.data, OrgDb = org.Mm.eg.db,
                                  type = "BP", organism = "mmu", # GO Biological Process
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, topn = 5) # Stricter cutoffs
      
      # Visualize
      pdf(file.path(RESULTS_DIR_SUBSET, 'clustergvis_heatmap_enrichment.pdf'), height = 10, width = 12)
      print(visCluster(object = st.data, plot.type = "both", # Heatmap and enrichment
                       markGenes = unique(top_markers_for_gvis$gene),
                       annoTerm.data = enrich.res,
                       sample.col = name2_cols, # Colors for 'name2' clusters
                       go.col = RColorBrewer::brewer.pal(5,"Set2")[1:5], # Colors for GO terms (example)
                       ctAnno.col = name2_cols, # Annotation bar colors
                       column_names_rot = 45))
      dev.off()
    } else {
      message("No significant top markers found for ClusterGVis after filtering.")
    }
  }, error = function(e) {
    message(paste("ClusterGVis failed:", e$message))
  })
}


# --- 8.8 Other Visualizations and Analyses from Original Script ---
message("--- Part 8.8: Additional Plots and Analyses ---")

# 8.8.1 GO Module Score Plot (if 'tnn' object and GO lists are defined)
# This part relied on an object 'tnn' and 'genelist' which are not clearly defined earlier in this segment.
go_term_file <- file.path(BASE_DIR, "../compare_development/GO/GO0002063_term_chondrocyte development_20250122_090327.txt")
if (file.exists(go_term_file) && "group_id" %in% colnames(related@meta.data)) {
  message("Plotting GO Module Scores...")
  tryCatch({
    matrixgenes_df <- read.table(go_term_file, header = TRUE, stringsAsFactors = FALSE) # Adjust if no header
    # Assuming second column contains gene symbols
    go_genelist <- list("ChondroDev" = unique(matrixgenes_df[,2])) # Name the list for AddModuleScore
    
    related <- AddModuleScore(object = related, features = go_genelist,
                              name = 'ChondroDevScore', assay = 'RNA', search=TRUE) # search=T if genes are not all present
    
    pdf(file.path(RESULTS_DIR_SUBSET, 'go_module_score_featuredimplot.pdf'), width = 8, height = 4)
    print(SCP::FeatureDimPlot(related, features = names(go_genelist)[1], # Plot by the score name (e.g., ChondroDevScore1)
                              pt.size = 1, ncol = 1, reduction = "umap.scvi_subset",
                              split.by = 'group_id', theme_light = FALSE))
    dev.off()
  }, error = function(e) {
    message(paste("GO Module Score plotting failed:", e$message))
  })
}


# --- 8.9 Final Cleanup ---
message("--- Part 8.9: Final Cleanup ---")
# rm(...) # Remove specific large intermediate objects if needed
gc()

sessionInfo()
message(paste0("Subset analysis script execution finished at ", Sys.time(), "! Results are in: ", RESULTS_DIR_SUBSET))