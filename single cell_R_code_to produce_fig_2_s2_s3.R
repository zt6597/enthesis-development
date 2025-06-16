# --------------------------------------------------------------------------
# Script: Developmental scRNA-seq Analysis (e.g., Tendon/Enthesis Development)
# Organization Date: 06-16-2025
# Description: This script processes entheses scRNA-seq data from different developmental
#              time points. It includes:
#              1. Loading data for early time points (e.g., e15, p07).
#              2. QC, DoubletFinder, SoupX correction.
#              3. Seurat SCT integration of these early time points.
#              4. Subsetting of relevant cell populations.
#              5. Re-integration of the subset.
#              6. Loading a reference dataset (e.g., p7-p28 annotated).
#              7. Label transfer from reference to the e15/p07 subset.
#              8. Merging the annotated e15/p07 data with the reference.
#              9. Full re-integration of the combined developmental timeline.
#              10. Clustering, visualization, and annotation of the final integrated object.
# --------------------------------------------------------------------------

# --- Environment Setup and Configuration ---
# Clear workspace (optional)
# rm(list=ls())
gc()

# Directories
BASE_DIR_DEV <- if (exists("project_base_dir_dev")) project_base_dir_dev else "~/compare_development/" # Main project dir
DATA_DIR_DEV <- file.path(BASE_DIR_DEV, "RC BTJ development data") # Input data dir
RESULTS_DIR_DEV <- file.path(BASE_DIR_DEV, "results_development_organized")
if (!dir.exists(RESULTS_DIR_DEV)) {
  dir.create(RESULTS_DIR_DEV, recursive = TRUE)
  message(paste("Created results directory:", RESULTS_DIR_DEV))
}
REFERENCE_RDS_PATH <- file.path(BASE_DIR_DEV, "rds file/normal_vs_unloaded_entheses_@_p7_p14_p28_annotated.rds") # Path to reference dataset

# File Patterns for Data Loading
PATTERN_10X_DEV <- "zt_iltered_feature_bc_matrix" # Original: 'zt_iltered_feature_bc_matrix' (typo?) -> corrected to 'filtered'
PATTERN_BGI_DEV <- "tp_filtered_feature_bc_matrix" # This was not used in the initial merge in the script

# QC Parameters (same as main script, adjust if different for developmental data)
MIN_FEATURE_RNA_QUANTILE_DEV <- 0.05
MAX_PERCENT_MT_DEV <- 10
MT_GENE_PATTERN_DEV <- "^mt-"
RIBO_GENE_PATTERN_DEV <- "^Rp[sl]"
HSP_GENE_PATTERN_DEV <- "^Hsp"
DISSOCIATION_GENES_DEV <- c(
  'Atf3','Zfp36','Fos','Dusp1','Egr1','Gadd45g','Jun','Junb','Jund',
  'Actg1','Ier2','Sat1','Ubc','Nr4a1','Egr2','Fosb','Ppp1r15a',
  'Mt1','Mt2','Errfi1','Ier3','Nfkbia','Zfp36l1','Cebpb','Cebpd','Btg2'
) %>% unique()

# DoubletFinder Parameters
DOUBLET_RATE_DEV <- 0.076
DF_ASSAY_DEV <- "SCT"
DF_PCS_DEV <- 1:50 # Original used 1:50
DF_RESOLUTION_DEV <- 0.4

# Seurat Integration Parameters
INTEGRATION_NFEATURES <- 4000
INTEGRATION_K_WEIGHT_INITIAL <- 100 # For e15_p07
INTEGRATION_K_WEIGHT_FULL <- 80    # For full timeline e15-4w
INTEGRATION_DIMS <- 1:50

# Clustering Parameters
INITIAL_CLUSTER_RES <- 0.4
SUBSET_CLUSTER_RES <- 0.2
FINAL_CLUSTER_RES <- 0.4

# File Names for RDS objects
BEFORE_QC_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_seurat_merged_before_QC.rds")
QC_LIST_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_seurat_list_after_QC.rds")
DOUBLETFINDER_LIST_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_seurat_list_after_doubletfinder.rds")
SOUPX_CORRECTED_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_e15_p07_after_soupX.rds")
INITIAL_INTEGRATED_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_e15_p07_SCT_integrated.rds")
SUBSET_INTEGRATED_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_e15_p07_subset_SCT_integrated.rds")
ANNOTATED_SUBSET_RDS_DEV <- file.path(RESULTS_DIR_DEV, "dev_e15_p07_subset_annotated.rds")
FINAL_INTEGRATED_TIMELINE_RDS <- file.path(RESULTS_DIR_DEV, "dev_e15_4w_enthesis_integrated_annotated.rds")

# Gene list files (developmental specific)
MT_GENES_FILE_DEV <- file.path(RESULTS_DIR_DEV, "cumulative_mitochondrial_genes_removed_dev.csv")
RIBO_GENES_FILE_DEV <- file.path(RESULTS_DIR_DEV, "cumulative_ribosomal_genes_removed_dev.csv")
HSP_GENES_FILE_DEV <- file.path(RESULTS_DIR_DEV, "cumulative_hsp_genes_removed_dev.csv")
DISS_GENES_FILE_DEV <- file.path(RESULTS_DIR_DEV, "cumulative_dissociation_genes_removed_dev.csv")
file.remove(MT_GENES_FILE_DEV, RIBO_GENES_FILE_DEV, HSP_GENES_FILE_DEV, DISS_GENES_FILE_DEV, showWarnings = FALSE)

# Options for `future` package to handle large objects
options(future.globals.maxSize = 4 * 1024^3) # 4 GiB, adjust as needed

# --- 9.1 Data Loading (Early Developmental Timepoints e.g., e15, p07) ---
message("--- Part 9.1: Loading Early Developmental Data ---")

# Define the helper functions `load_scrna_data` and `filter_and_record_genes`
# Helper function placeholder (copy from previous script if not sourced)
load_scrna_data <- function(data_dir, folder_pattern, data_type_prefix, gene_col_for_bgi = 2) { # gene.column = 2 for 10x default
  folders <- list.files(data_dir, pattern = folder_pattern, full.names = FALSE)
  if (length(folders) == 0) {
    warning(paste("No folders found matching pattern:", folder_pattern, "in directory:", data_dir))
    return(list())
  }
  sce_list <- lapply(folders, function(folder_name) {
    message(paste("Processing folder:", file.path(data_dir, folder_name)))
    counts_path <- file.path(data_dir, folder_name)
    current_counts <- Read10X(counts_path, gene.column = if (grepl("bgi", folder_pattern, ignore.case = TRUE)) gene_col_for_bgi else 2)
    seurat_obj <- CreateSeuratObject(counts = current_counts, project = folder_name)
    return(seurat_obj)
  })
  names(sce_list) <- folders
  return(sce_list)
}

filter_and_record_genes <- function(seurat_obj, gene_pattern = NULL, gene_list = NULL, output_file, gene_set_name) {
  rna_assay <- seurat_obj@assays$RNA
  if (is.null(rna_assay)) {
    message(sprintf("RNA assay not found for sample %s. Skipping %s gene filtering.", unique(seurat_obj$orig.ident), gene_set_name))
    return(seurat_obj)
  }
  all_genes <- rownames(rna_assay)
  genes_to_remove_indices <- c()
  if (!is.null(gene_pattern)) genes_to_remove_indices <- c(genes_to_remove_indices, grep(gene_pattern, all_genes, ignore.case = TRUE))
  if (!is.null(gene_list)) genes_to_remove_indices <- c(genes_to_remove_indices, which(all_genes %in% gene_list))
  genes_to_remove_indices <- unique(genes_to_remove_indices)
  if (length(genes_to_remove_indices) > 0) {
    genes_identified <- all_genes[genes_to_remove_indices]
    if (length(genes_identified) > 0) {
      if (!file.exists(output_file)) write.table(data.frame(Gene = character()), output_file, sep = ",", col.names = TRUE, row.names = FALSE, quote=FALSE)
      existing_genes <- tryCatch(read.csv(output_file, stringsAsFactors = FALSE)$Gene, error = function(e) character(0))
      new_genes_to_write <- setdiff(genes_identified, existing_genes)
      if(length(new_genes_to_write) > 0) write.table(data.frame(Gene = new_genes_to_write), output_file, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE, quote=FALSE)
    }
    features_to_keep <- setdiff(all_genes, genes_identified)
    if (length(features_to_keep) < length(all_genes)) {
      seurat_obj@assays$RNA <- subset(rna_assay, features = features_to_keep)
      message(sprintf("Removed %d %s genes from RNA assay of sample %s. %d features remaining.", length(genes_identified), gene_set_name, unique(seurat_obj$orig.ident), length(features_to_keep)))
    }
  } else { message(sprintf("No %s genes to remove found in sample %s.", gene_set_name, unique(seurat_obj$orig.ident))) }
  return(seurat_obj)
}

# Load "10X-like" data (zt_filtered)
sce_list_dev_10x <- load_scrna_data(DATA_DIR_DEV, folder_pattern = PATTERN_10X_DEV, data_type_prefix = "10X")

# Load "BGI-like" data (tp_filtered) - Note: original script did not merge this list initially.
# sce_list_dev_bgi <- load_scrna_data(DATA_DIR_DEV, folder_pattern = PATTERN_BGI_DEV, data_type_prefix = "BGI")

# Initial merge (original script only merged sceList10x items)
if (length(sce_list_dev_10x) == 0) {
  stop("No developmental data loaded (sce_list_dev_10x is empty). Check DATA_DIR_DEV and PATTERN_10X_DEV.")
}
if (length(sce_list_dev_10x) == 1) {
  seurat_merged_raw_dev <- sce_list_dev_10x[[1]]
} else {
  seurat_merged_raw_dev <- merge(sce_list_dev_10x[[1]],
                                 y = sce_list_dev_10x[2:length(sce_list_dev_10x)],
                                 add.cell.ids = names(sce_list_dev_10x))
}
seurat_merged_raw_dev@project.name <- "Compare_Development_Raw"
saveRDS(seurat_merged_raw_dev, BEFORE_QC_RDS_DEV)

# --- 9.2 Quality Control (per sample) ---
message("--- Part 9.2: Quality Control for Developmental Samples ---")
# The original script's QC loop was on `sceList10x`.
sce_list_filtered_dev <- lapply(X = sce_list_dev_10x, FUN = function(x) { # Using sce_list_dev_10x as per original
  current_ident <- unique(x$orig.ident)
  message(paste("Performing QC for sample:", current_ident))
  
  # Metadata extraction (adjust substr indices based on your folder names for dev data)
  # Example: "e15_myrun_zt_filtered..." or "p07_another_zt_filtered..."
  x$sample_full_id <- current_ident # Keep full original ident
  x$time_point <- substr(current_ident, 1, 3) # e.g., "e15", "p07"
  # Add other metadata like 'group' or 'tech' if applicable and derivable from folder name
  # x$group_dev  <- substr(current_ident, ?, ?)
  
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = MT_GENE_PATTERN_DEV, assay = "RNA")
  x[["percent.ribo"]] <- PercentageFeatureSet(x, pattern = RIBO_GENE_PATTERN_DEV, assay = "RNA")
  
  nfeature_rna_low_threshold <- quantile(x$nFeature_RNA, probs = MIN_FEATURE_RNA_QUANTILE_DEV, na.rm = TRUE)
  cells_before_filter <- ncol(x)
  x <- subset(x, subset = nFeature_RNA > nfeature_rna_low_threshold & percent.mt < MAX_PERCENT_MT_DEV)
  message(sprintf("Filtered cells for %s: %d to %d (nFeature_RNA > %.0f, percent.mt < %d)",
                  current_ident, cells_before_filter, ncol(x), nfeature_rna_low_threshold, MAX_PERCENT_MT_DEV))
  if (ncol(x) == 0) return(NULL)
  
  x <- filter_and_record_genes(x, gene_pattern = MT_GENE_PATTERN_DEV, output_file = MT_GENES_FILE_DEV, gene_set_name = "Mitochondrial")
  x <- filter_and_record_genes(x, gene_pattern = RIBO_GENE_PATTERN_DEV, output_file = RIBO_GENES_FILE_DEV, gene_set_name = "Ribosomal")
  x <- filter_and_record_genes(x, gene_pattern = HSP_GENE_PATTERN_DEV, output_file = HSP_GENES_FILE_DEV, gene_set_name = "Heat Shock Protein")
  x <- filter_and_record_genes(x, gene_list = DISSOCIATION_GENES_DEV, output_file = DISS_GENES_FILE_DEV, gene_set_name = "Dissociation-induced")
  
  return(x)
})
sce_list_filtered_dev <- sce_list_filtered_dev[!sapply(sce_list_filtered_dev, is.null)]
if(length(sce_list_filtered_dev) == 0) stop("All dev samples filtered out during QC.")
saveRDS(sce_list_filtered_dev, QC_LIST_RDS_DEV)
#rm(sceList10x); gc() # sce_list_dev_10x is used by before_QC for SoupX

# --- 9.3 Doublet Detection and Removal (per sample) ---
message("--- Part 9.3: Doublet Detection for Developmental Samples ---")
# sce_list_filtered_dev <- readRDS(QC_LIST_RDS_DEV) # Optional load

sce_list_after_doubletfinder_dev <- lapply(X = sce_list_filtered_dev, FUN = function(x) {
  current_ident <- unique(x$orig.ident)
  message(paste("Running DoubletFinder for dev sample:", current_ident))
  x <- SCTransform(x, assay="RNA", new.assay.name = "SCT", vst.flavor = "v2", variable.features.n = 4000, verbose = FALSE)
  x <- RunPCA(x, assay = DF_ASSAY_DEV, npcs = max(DF_PCS_DEV), verbose = FALSE)
  x <- RunUMAP(x, assay = DF_ASSAY_DEV, dims = DF_PCS_DEV, reduction.name = "umap.DF", reduction.key = "UMAPDF_", verbose = FALSE)
  x <- FindNeighbors(x, assay = DF_ASSAY_DEV, dims = DF_PCS_DEV, verbose = FALSE)
  x <- FindClusters(x, resolution = DF_RESOLUTION_DEV, verbose = FALSE)
  
  sweep.res.list <- paramSweep_v3(x, PCs = DF_PCS_DEV, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  optimal_pK <- if (nrow(bcmvn) > 0 && any(!is.na(bcmvn$BCmetric))) as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])) else 0.09
  
  annotations <- x$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(DOUBLET_RATE_DEV * ncol(x)) 
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  x <- doubletFinder_v3(x, PCs = DF_PCS_DEV, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  df_classification_col <- colnames(x@meta.data)[grep(paste0("^DF.classifications_0.25_", optimal_pK), colnames(x@meta.data))]
  if (length(df_classification_col)==0) df_classification_col <- colnames(x@meta.data)[grep("^DF.classifications", colnames(x@meta.data))][1] # Fallback
  
  x$doublet_classification <- x@meta.data[[df_classification_col]]
  singlets <- subset(x, subset = doublet_classification == "Singlet")
  message(sprintf("Identified doublets for %s. Kept %d singlets out of %d cells.", current_ident, ncol(singlets), ncol(x)))
  
  # Return object with RNA assay of singlets (original counts)
  singlet_obj_rna_only <- CreateSeuratObject(counts = singlets@assays$RNA@counts, meta.data = singlets@meta.data, project = unique(singlets$orig.ident))
  return(singlet_obj_rna_only)
})
sce_list_after_doubletfinder_dev <- sce_list_after_doubletfinder_dev[!sapply(sce_list_after_doubletfinder_dev, is.null)]
if(length(sce_list_after_doubletfinder_dev) == 0) stop("All dev samples filtered out after DoubletFinder.")
saveRDS(sce_list_after_doubletfinder_dev, DOUBLETFINDER_LIST_RDS_DEV)
rm(sce_list_filtered_dev); gc()

# --- 9.4 Ambient RNA Correction with SoupX ---
message("--- Part 9.4: SoupX Correction for Developmental Samples ---")
# sce_list_after_doubletfinder_dev <- readRDS(DOUBLETFINDER_LIST_RDS_DEV) # Optional load

# Manual metadata correction (as in original script)
# This assumes sce_list_after_doubletfinder_dev[[2]] is p07, and sce_list_dev_10x[[2]] is p07
if (length(sce_list_after_doubletfinder_dev) >= 2) {
  sce_list_after_doubletfinder_dev[[2]]$time_point <- 'p07'
}
if (length(sce_list_dev_10x) >= 2) {
  sce_list_dev_10x[[2]]$time_point <- 'p07' 
}

# Merge for toc (after QC, after DoubletFinder)
if (length(sce_list_after_doubletfinder_dev) == 1) {
  seurat_for_toc_dev <- sce_list_after_doubletfinder_dev[[1]]
} else {
  sample_names_toc <- sapply(sce_list_after_doubletfinder_dev, function(s) unique(s$orig.ident))
  seurat_for_toc_dev <- merge(sce_list_after_doubletfinder_dev[[1]],
                              y = sce_list_after_doubletfinder_dev[2:length(sce_list_after_doubletfinder_dev)],
                              add.cell.ids = make.unique(sample_names_toc))
}
toc_dev <- seurat_for_toc_dev@assays$RNA@counts

# Merge for tod (before QC, raw counts - using sce_list_dev_10x as per original script's logic)
if (length(sce_list_dev_10x) == 1) {
  seurat_for_tod_dev <- sce_list_dev_10x[[1]]
} else {
  sample_names_tod <- sapply(sce_list_dev_10x, function(s) unique(s$orig.ident))
  seurat_for_tod_dev <- merge(sce_list_dev_10x[[1]],
                              y = sce_list_dev_10x[2:length(sce_list_dev_10x)],
                              add.cell.ids = make.unique(sample_names_tod))
}
tod_dev <- seurat_for_tod_dev@assays$RNA@counts

# Align genes and cells for SoupX
common_genes_dev <- intersect(rownames(tod_dev), rownames(toc_dev))
tod_dev <- tod_dev[common_genes_dev, ]
toc_dev <- toc_dev[common_genes_dev, ]
common_cells_dev <- intersect(colnames(tod_dev), colnames(toc_dev))
toc_dev <- toc_dev[, common_cells_dev]
seurat_for_toc_dev <- subset(seurat_for_toc_dev, cells = common_cells_dev)

# Quick clustering for SoupX
quick_cluster_obj_dev <- seurat_for_toc_dev
DefaultAssay(quick_cluster_obj_dev) <- "RNA"
quick_cluster_obj_dev <- NormalizeData(quick_cluster_obj_dev, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
  RunPCA(features = VariableFeatures(.), npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:15, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
soup_clusters_dev <- setNames(quick_cluster_obj_dev$seurat_clusters, colnames(quick_cluster_obj_dev))

# Run SoupX
sc_dev <- SoupChannel(tod_dev, toc_dev, calcSoupProfile = TRUE)
sc_dev <- setClusters(sc_dev, soup_clusters_dev)
sc_dev <- autoEstCont(sc_dev, verbose = TRUE)
out_counts_dev <- adjustCounts(sc_dev, roundToInt = TRUE)

# Create Seurat object with corrected counts
seurat_corrected_dev <- CreateSeuratObject(counts = out_counts_dev, meta.data = seurat_for_toc_dev@meta.data)
seurat_corrected_dev@project.name <- "Dev_e15_p07_SoupX"
# Factor time_point (e.g. e15, p07)
seurat_corrected_dev$time_point <- factor(seurat_corrected_dev$time_point, levels = intersect(c('e15','p07'), unique(seurat_corrected_dev$time_point)))
saveRDS(seurat_corrected_dev, SOUPX_CORRECTED_RDS_DEV)
rm(seurat_for_tod_dev, seurat_for_toc_dev, tod_dev, toc_dev, quick_cluster_obj_dev, sc_dev, out_counts_dev, sce_list_after_doubletfinder_dev, sce_list_dev_10x); gc()

# --- 9.5 Seurat Integration (SCTransform CCA) of Early Timepoints (e.g., e15, p07) ---
message("--- Part 9.5: Seurat Integration of Early Developmental Timepoints ---")
seurat_corrected_dev <- readRDS(SOUPX_CORRECTED_RDS_DEV) # Load SoupX corrected data

dev_list <- SplitObject(seurat_corrected_dev, split.by = "time_point")
dev_list <- lapply(X = dev_list, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", vars.to.regress = "percent.mt", variable.features.n = INTEGRATION_NFEATURES, verbose = FALSE)
  return(x)
})
features_dev <- SelectIntegrationFeatures(object.list = dev_list, nfeatures = INTEGRATION_NFEATURES)
dev_list <- PrepSCTIntegration(object.list = dev_list, anchor.features = features_dev)
integrate_anchors_dev <- FindIntegrationAnchors(object.list = dev_list, normalization.method = "SCT", anchor.features = features_dev)
seurat_integrated_early <- IntegrateData(anchorset = integrate_anchors_dev, normalization.method = "SCT", k.weight = INTEGRATION_K_WEIGHT_INITIAL)

DefaultAssay(seurat_integrated_early) <- 'integrated'
seurat_integrated_early <- RunPCA(seurat_integrated_early, verbose = FALSE)
seurat_integrated_early <- FindNeighbors(seurat_integrated_early, dims = INTEGRATION_DIMS)
seurat_integrated_early <- FindClusters(seurat_integrated_early, resolution = INITIAL_CLUSTER_RES)
seurat_integrated_early <- RunUMAP(seurat_integrated_early, dims = INTEGRATION_DIMS, reduction.name = "umap_early")

DimPlot(seurat_integrated_early, reduction = "umap_early", group.by = "seurat_clusters", label = TRUE) + ggtitle("e15_p07 Integrated UMAP")
saveRDS(seurat_integrated_early, INITIAL_INTEGRATED_RDS_DEV)
# seurat_integrated_early <- readRDS(INITIAL_INTEGRATED_RDS_DEV)

# --- 9.6 Subset Relevant Cells and Re-integrate ---
message("--- Part 9.6: Subsetting and Re-integrating Relevant Cell Population ---")
Idents(seurat_integrated_early) <- 'seurat_clusters'
clusters_to_subset_early <- c('1','2','6','9')
related_early_subset <- subset(seurat_integrated_early, idents = clusters_to_subset_early)

# Re-create with RNA counts for new integration split by "sample_full_id" (original orig.ident)
related_early_diet <- CreateSeuratObject(counts = related_early_subset@assays$RNA@counts, meta.data = related_early_subset@meta.data)
related_early_list <- SplitObject(related_early_diet, split.by = "sample_full_id") 

related_early_list <- lapply(X = related_early_list, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", vars.to.regress = "percent.mt", variable.features.n = INTEGRATION_NFEATURES, verbose = FALSE)
  return(x)
})
features_related_early <- SelectIntegrationFeatures(object.list = related_early_list, nfeatures = INTEGRATION_NFEATURES)
related_early_list <- PrepSCTIntegration(object.list = related_early_list, anchor.features = features_related_early)
integrate_anchors_related_early <- FindIntegrationAnchors(object.list = related_early_list, normalization.method = "SCT", anchor.features = features_related_early)
related_early_integrated <- IntegrateData(anchorset = integrate_anchors_related_early, normalization.method = "SCT", k.weight = INTEGRATION_K_WEIGHT_INITIAL)

DefaultAssay(related_early_integrated) <- 'integrated'
related_early_integrated <- RunPCA(related_early_integrated, verbose = FALSE)
related_early_integrated <- FindNeighbors(related_early_integrated, dims = INTEGRATION_DIMS)
related_early_integrated <- FindClusters(related_early_integrated, resolution = SUBSET_CLUSTER_RES)
related_early_integrated <- RunUMAP(related_early_integrated, dims = INTEGRATION_DIMS, reduction.name = "umap_related_early")

Idents(related_early_integrated) <- "seurat_clusters"
clusters_to_remove_subset <- c('4','5') # Adjust as needed
related_early_integrated <- subset(related_early_integrated, idents = clusters_to_remove_subset, invert = TRUE)

DimPlot(related_early_integrated, reduction = "umap_related_early", group.by = "seurat_clusters", label = TRUE) + ggtitle("e15_p07 Subset Re-integrated UMAP")
saveRDS(related_early_integrated, SUBSET_INTEGRATED_RDS_DEV)
# related_early_integrated <- readRDS(SUBSET_INTEGRATED_RDS_DEV)

# --- 9.7 Reference Loading and Label Transfer ---
message("--- Part 9.7: Reference Loading and Label Transfer ---")
ref_later_timeline <- readRDS(REFERENCE_RDS_PATH)
# Subset reference to 'nom' group as in original script
Idents(ref_later_timeline) <- 'name1' 
if ("nom" %in% levels(Idents(ref_later_timeline))) {
  ref_later_timeline <- subset(ref_later_timeline, idents = "nom")
} else {
  warning("'nom' group not found in reference 'name1' idents. Using full reference.")
}

DefaultAssay(ref_later_timeline) <- "RNA"
ref_later_timeline <- SCTransform(ref_later_timeline, assay = "RNA", new.assay.name = "SCT_ref",
                                  vst.flavor = "v2", vars.to.regress = "percent.mt",
                                  variable.features.n = INTEGRATION_NFEATURES, verbose = FALSE)
# Reference needs PCA on its SCT assay
DefaultAssay(ref_later_timeline) <- "SCT_ref"
ref_later_timeline <- RunPCA(ref_later_timeline, assay = "SCT_ref", reduction.name = "pca_ref", verbose = FALSE)


# Prepare query (related_early_integrated)
DefaultAssay(related_early_integrated) <- "RNA" # Anchors work on RNA, normalization.method="SCT" handles it

related_early_integrated <- SCTransform(related_early_integrated, assay="RNA", new.assay.name="SCT_query",
                                        vst.flavor = "v2", vars.to.regress = "percent.mt",
                                        variable.features.n = INTEGRATION_NFEATURES, verbose = FALSE)
DefaultAssay(related_early_integrated) <- "SCT_query" # Ensure this is used by query

transfer_anchors <- FindTransferAnchors(reference = ref_later_timeline, query = related_early_integrated,
                                        reference.assay = "SCT_ref", query.assay = "SCT_query", # Specify SCT assays
                                        normalization.method = "SCT", # Already SCTed, but good practice
                                        dims = INTEGRATION_DIMS, reference.reduction = "pca_ref")

predictions_dev <- TransferData(anchorset = transfer_anchors, refdata = ref_later_timeline$name2, # Transfer 'name2' labels from reference
                                dims = INTEGRATION_DIMS) # Use the same dims
related_early_annotated <- AddMetaData(related_early_integrated, metadata = predictions_dev)
related_early_annotated$name2 <- related_early_annotated$predicted.id # Store predicted labels in 'name2'

saveRDS(related_early_annotated, ANNOTATED_SUBSET_RDS_DEV)
# related_early_annotated <- readRDS(ANNOTATED_SUBSET_RDS_DEV)

# --- 9.8 Merge Annotated Early Data with Reference for Full Timeline ---
message("--- Part 9.8: Merging Annotated Early Data with Reference ---")
# Ensure 'name2' exists in ref_later_timeline if it's the common annotation column
if (!"name2" %in% colnames(ref_later_timeline@meta.data)) {
  warning("'name2' column not found in reference. Label transfer might not have completed as expected or ref needs 'name2'.")
}
# Ensure 'time_point' is consistent, reference might have 'time'
if ("time" %in% colnames(ref_later_timeline@meta.data) && !"time_point" %in% colnames(ref_later_timeline@meta.data)) {
  ref_later_timeline$time_point <- ref_later_timeline$time
}

# Select common metadata columns to avoid issues, ensure 'RNA' assay is present with counts
common_meta_cols <- intersect(colnames(ref_later_timeline@meta.data), colnames(related_early_annotated@meta.data))
ref_for_merge <- subset(ref_later_timeline, features = rownames(ref_later_timeline[["RNA"]])) # Ensure RNA assay is primary
ref_for_merge@meta.data <- ref_for_merge@meta.data[, common_meta_cols, drop=FALSE]

query_for_merge <- subset(related_early_annotated, features = rownames(related_early_annotated[["RNA"]]))
query_for_merge@meta.data <- query_for_merge@meta.data[, common_meta_cols, drop=FALSE]

# Check for unique cell names before merge
if (any(colnames(ref_for_merge) %in% colnames(query_for_merge))) {
  message("Cell names overlap, adding prefixes for merge.")
  ref_for_merge <- RenameCells(ref_for_merge, add.cell.id = "ref")
  query_for_merge <- RenameCells(query_for_merge, add.cell.id = "query")
}


# Create a diet object for merging based on RNA counts and essential metadata
meta_ref <- ref_later_timeline@meta.data
counts_ref <- ref_later_timeline@assays$RNA@counts
ref_diet_for_merge <- CreateSeuratObject(counts_ref, meta.data = meta_ref)

meta_query <- related_early_annotated@meta.data
counts_query <- related_early_annotated@assays$RNA@counts
query_diet_for_merge <- CreateSeuratObject(counts_query, meta.data = meta_query)

# Add unique cell IDs if not already done
ref_diet_for_merge <- RenameCells(ref_diet_for_merge, new.names = paste0("ref_", colnames(ref_diet_for_merge)))
query_diet_for_merge <- RenameCells(query_diet_for_merge, new.names = paste0("query_", colnames(query_diet_for_merge)))

full_timeline_raw_merged <- merge(ref_diet_for_merge, y = query_diet_for_merge)
full_timeline_raw_merged@project.name <- "Full_Dev_Timeline_PreIntegration"

# Integrate the full timeline (e.g. e15, p07, 1w, 2w, 4w)
# Split by original sample identifier for robust integration
# `sample_full_id` was created during QC for dev data. Reference needs a similar unique sample ID column.
# Assuming reference has `sample_id` or `orig.ident` that can serve this purpose.
if (!"sample_full_id" %in% colnames(full_timeline_raw_merged@meta.data)) {
  full_timeline_raw_merged$sample_full_id <- full_timeline_raw_merged$orig.ident
}

full_timeline_list <- SplitObject(full_timeline_raw_merged, split.by = "sample_full_id")

full_timeline_list <- lapply(X = full_timeline_list, FUN = function(x) {
  # Ensure percent.mt is calculated if not present for all objects in list
  if (!"percent.mt" %in% colnames(x@meta.data)) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = MT_GENE_PATTERN_DEV, assay = "RNA")
  }
  x <- SCTransform(x, assay="RNA", new.assay.name="SCT", vst.flavor = "v2",
                   vars.to.regress = "percent.mt", variable.features.n = INTEGRATION_NFEATURES, verbose = FALSE)
  return(x)
})
features_full_timeline <- SelectIntegrationFeatures(object.list = full_timeline_list, nfeatures = INTEGRATION_NFEATURES)
full_timeline_list <- PrepSCTIntegration(object.list = full_timeline_list, anchor.features = features_full_timeline)
integrate_anchors_full_timeline <- FindIntegrationAnchors(object.list = full_timeline_list, normalization.method = "SCT", anchor.features = features_full_timeline)
enthesis_timeline_integrated <- IntegrateData(anchorset = integrate_anchors_full_timeline, normalization.method = "SCT", k.weight = INTEGRATION_K_WEIGHT_FULL)

DefaultAssay(enthesis_timeline_integrated) <- 'integrated'
enthesis_timeline_integrated <- RunPCA(enthesis_timeline_integrated, verbose = FALSE)
enthesis_timeline_integrated <- FindNeighbors(enthesis_timeline_integrated, dims = INTEGRATION_DIMS)
enthesis_timeline_integrated <- FindClusters(enthesis_timeline_integrated, resolution = FINAL_CLUSTER_RES)
enthesis_timeline_integrated <- RunUMAP(enthesis_timeline_integrated, dims = INTEGRATION_DIMS, reduction.name = "umap_full_timeline")

DimPlot(enthesis_timeline_integrated, reduction = "umap_full_timeline", group.by = "time_point", label = TRUE) + ggtitle("Full Developmental Timeline (e15-4w) UMAP")

# --- 9.9 Final Annotation and Visualization of Full Timeline ---
message("--- Part 9.9: Final Annotation and Visualization of Full Timeline ---")
# The object 'enthesis_timeline_integrated' now contains cells from e15, p07 (query), and p7, 1w, 2w, 4w (reference).
# 'name2' was transferred to e15/p07 and should exist for reference cells.

# Factor 'time_point'
expected_time_levels <- c('e15','p07','1w','2w','4w') # p01 was used in original script, mapping from p07?
enthesis_timeline_integrated$time_point <- factor(enthesis_timeline_integrated$time_point, levels = intersect(expected_time_levels, unique(enthesis_timeline_integrated$time_point)))

# Use existing 'name2' for plotting, ensure it's factored
name2_levels_final <- unique(enthesis_timeline_integrated$name2)
name2_levels_final <- name2_levels_final[!is.na(name2_levels_final)]
enthesis_timeline_integrated$name2 <- factor(enthesis_timeline_integrated$name2, levels = intersect(c('Enthesis chondrocytes','Tenocytes','Enthesis progenitors','Tenoblasts','Mesenchymal progenitors', 'related chondrocytes', 'related progenitors'), name2_levels_final)) # Add levels from prediction

name2_cols_final <- ggsci::pal_locuszoom("default")(length(levels(enthesis_timeline_integrated$name2)))
names(name2_cols_final) <- levels(enthesis_timeline_integrated$name2)

DimPlot(enthesis_timeline_integrated, reduction = "umap_full_timeline", group.by = "name2", label = TRUE, cols = name2_cols_final, pt.size=0.5) + ggtitle("Full Timeline by 'name2' Annotation")

# Scale RNA data for feature plots and other RNA-based analyses
DefaultAssay(enthesis_timeline_integrated) <- 'RNA'
enthesis_timeline_integrated <- NormalizeData(enthesis_timeline_integrated, verbose = FALSE)
enthesis_timeline_integrated <- ScaleData(enthesis_timeline_integrated, features = rownames(enthesis_timeline_integrated),
                                          vars.to.regress = c('percent.mt','nCount_RNA','nFeature_RNA'), verbose = TRUE) # As in original

saveRDS(enthesis_timeline_integrated, FINAL_INTEGRATED_TIMELINE_RDS)
# enthesis_timeline_integrated <- readRDS(FINAL_INTEGRATED_TIMELINE_RDS)

# FeaturePlots
genes_to_plot_final <- c('Scx','Sox9','Col2a1','Tnn','Tnmd','Ly6a','Clec3a','Runx2')
FeaturePlot(enthesis_timeline_integrated, features = genes_to_plot_final, reduction = "umap_full_timeline", order=TRUE, pt.size=0.3)

# Violin Plot
VlnPlot(enthesis_timeline_integrated, features = 'Tnn', assay = 'RNA', group.by = 'time_point', pt.size=0)

# CellStatPlot
if (requireNamespace("scplotter", quietly = TRUE) && "name2" %in% colnames(enthesis_timeline_integrated@meta.data)) {
  pdf(file.path(RESULTS_DIR_DEV, "cellstat_full_timeline_name2_by_time.pdf"), width=8, height=6)
  print(scplotter::CellStatPlot(enthesis_timeline_integrated, group_by = "name2", stat_by = "time_point",
                                palcolor = name2_cols_final, position = "stack", aspect.ratio = 1))
  dev.off()
}

# --- 9.10 Final Cleanup ---
message("--- Part 9.10: Final Cleanup ---")
gc()

sessionInfo()
message(paste0("Developmental analysis script execution finished at ", Sys.time(), "! Results are in: ", RESULTS_DIR_DEV))

# Description: This script performs advanced analyses on the fully integrated
#              developmental Seurat object ('enthesis'). It includes:
#              1. Differential Gene Expression (DEG) analysis for tendon and enthesis subclusters.
#              2. Trajectory inference using Slingshot (via SCP) on Diffusion Map embedding.
#              3. Dynamic gene feature analysis along trajectories.
#              4. GO enrichment analysis of DEGs.
#              5. Subsetting 'chondrogenesis' related cells for focused analysis.
#              6. Trajectory analysis with GeneTrajectory package.
#              7. Data export to H5Seurat and AnnData formats.
#              8. CytoTRACE2 analysis (potentially on a different object 'tnn').
# Prerequisite: An integrated Seurat object named 'enthesis' from the previous
#               developmental script (e.g., from 'e15_4w_enthesis_SCT_RNA_scalded.rds').
# --------------------------------------------------------------------------
# Ensure 'enthesis' object is loaded
if (!exists("enthesis") || !inherits(enthesis, "Seurat")) {
  warning("'enthesis' Seurat object not found. Please load it first, e.g., from 'e15_4w_enthesis_SCT_RNA_scalded.rds'")
  # Example: enthesis <- readRDS(file.path(RESULTS_DIR_DEV, "e15_4w_enthesis_SCT_RNA_scalded.rds"))
  # This script will fail if 'enthesis' is not loaded.
}

# Parameters
DEG_LOG2FC_THRESHOLD <- 1
DEG_TOP_N_PER_GROUP <- 5
TRAJECTORY_START_NODE <- 'Mesenchymal progenitors'
TRAJECTORY_END_NODES <- c('Enthesis chondrocytes','Tenocytes')
TRAJECTORY_REDUCTION <- "dm" # Diffusion Map for Slingshot
DYN_FEAT_N_CANDIDATES <- 200
GENETRAJ_VAR_FEATURES <- 1000
GENETRAJ_EXPR_PERCENT_MIN <- 0.2
GENETRAJ_EXPR_PERCENT_MAX <- 0.5
GENETRAJ_META_CELLS_N <- 500
GENETRAJ_K <- 10 # For GetGraphDistance and GetGeneEmbedding
GENETRAJ_TRAJECTORIES_N <- 3

# File Names
ENTHESIS_DEGS_PLOT_PDF <- file.path(RESULTS_DIR_DEV, "enthesis_top_DEGs_heatmap.pdf")
ENTHESIS_LINEAGE_PLOT_PDF <- file.path(RESULTS_DIR_DEV, "enthesis_name2_lineage_dm.pdf")
ENTHESIS_DYNFEAT_PLOT_PDF <- file.path(RESULTS_DIR_DEV, "enthesis_name2_lineage_dynamic_genes.pdf")
ENTHESIS_DYNHEATMAP_GO_PDF <- file.path(RESULTS_DIR_DEV, "enthesis_name2_dynamic_heatmap_GO.pdf")
CHONDROGENESIS_SUBSET_RDS <- file.path(RESULTS_DIR_DEV, "chondrogenesis_subset_scaled.rds")
GENETRAJ_SCATTER3D_PDF <- file.path(RESULTS_DIR_DEV, "genetraj_chondro_scatter3d.pdf")
GENETRAJ_FEATUREPLOT_PDF <- file.path(RESULTS_DIR_DEV, "genetraj_chondro_bin_featureplots.pdf")
ENTHESIS_H5SEURAT_OUT <- file.path(RESULTS_DIR_DEV, "enthesis_final.h5Seurat")
ENTHESIS_H5AD_OUT_RNA <- file.path(RESULTS_DIR_DEV, "enthesis_final_RNA.h5ad") # For RNA assay data
ENTHESIS_H5AD_OUT_SCVI <- file.path(BASE_DIR_DEV, "scvi/enthesis_name2_refference.h5ad") # Original path for sceasy
CYTOTRACE_PLOTS_PDF <- file.path(RESULTS_DIR_DEV, "cytotrace_tnn_plots.pdf") # Assuming 'tnn' object

# Color palettes (ensure name2_cols is defined from previous script if needed)
if (!exists("name2_cols") && "name2" %in% colnames(enthesis@meta.data)) {
  name2_levels_final <- levels(enthesis$name2)
  name2_cols <- ggsci::pal_locuszoom("default")(length(name2_levels_final))
  names(name2_cols) <- name2_levels_final
}
if (!exists("sample_cols") && "sample_id" %in% colnames(enthesis@meta.data)) { # Example
  sample_cols <- RColorBrewer::brewer.pal(n=min(8, length(unique(enthesis$sample_id))), name="Set2")
  names(sample_cols) <- unique(enthesis$sample_id)
}


# Python Environment for GeneTrajectory (adjust if different)
GENETRAJ_PYTHON_PATH <- "/home/whc/miniconda3/envs/genetrajectory/bin/python"
GENETRAJ_VENV_PATH <- "/home/whc/miniconda3/envs/genetrajectory"

# --- 10.1 Differential Gene Expression (DEG) for 'name2' cell types ---
message("--- Part 10.1: DEG Analysis for 'name2' cell types ---")
if ("name2" %in% colnames(enthesis@meta.data)) {
  # Assuming DEGs were already computed and stored in @tools by SCP::RunDEtest from previous script.
  # If not, run it here:
  if (!("DEtest_name2" %in% names(enthesis@tools)) || !("AllMarkers_wilcox" %in% names(enthesis@tools$DEtest_name2))) {
    message("RunDEtest results for 'name2' not found. Running SCP::RunDEtest...")
    # Ensure correct assay is active for DE testing (e.g., SCT or integrated)
    # DefaultAssay(enthesis) <- "SCT" # Or 'integrated' if appropriate
    # enthesis <- SCP::RunDEtest(srt = enthesis, group_by = "name2", assay = DefaultAssay(enthesis),
    #                            fc.threshold = log2(1.2), only.pos = FALSE) # Example thresholds
  }
  
  if ("DEtest_name2" %in% names(enthesis@tools) && "AllMarkers_wilcox" %in% names(enthesis@tools$DEtest_name2)) {
    degs_name2 <- enthesis@tools$DEtest_name2$AllMarkers_wilcox
    
    top_degs_name2 <- degs_name2 %>%
      dplyr::filter(avg_log2FC > DEG_LOG2FC_THRESHOLD, p_val_adj < 0.05) %>% # Add p-value filter
      dplyr::group_by(group1) %>% # group1 is the cluster name
      dplyr::slice_max(n = DEG_TOP_N_PER_GROUP, order_by = avg_log2FC) %>%
      dplyr::ungroup()
    
    if (nrow(top_degs_name2) > 0) {
      features_for_heatmap <- unique(top_degs_name2$gene)
      pdf(ENTHESIS_DEGS_PLOT_PDF, height = max(4, 0.2*length(features_for_heatmap)), width = 8)
      # Using FeatureStatPlot as dotplot/heatmap combo
      print(scplotter::FeatureStatPlot(enthesis, features = features_for_heatmap, group_by = "name2", 
                                       plot_type = "dotheatmap", 
                                       palcolor = name2_cols, 
                                       colorscale_limits = c(-2,2), 
                                       add_bg = FALSE, flip = TRUE,
                                       show_row_names = TRUE, cluster_rows = FALSE))
      dev.off()
    } else {
      message("No DEGs found meeting the criteria for heatmap plotting.")
    }
  } else {
    message("DEG results for 'name2' not found in enthesis@tools. Skipping DEG heatmap.")
  }
} else {
  message("Metadata column 'name2' not found. Skipping DEG analysis.")
}

# --- 10.2 Trajectory Analysis with Slingshot (SCP) ---
message("--- Part 10.2: Trajectory Analysis with Slingshot (SCP) ---")
if ("name2" %in% colnames(enthesis@meta.data) && TRAJECTORY_REDUCTION %in% Reductions(enthesis)) {
  options(timeout=600000000000000000000) # Extremely large timeout
  
  if (TRAJECTORY_START_NODE %in% levels(enthesis$name2) && all(TRAJECTORY_END_NODES %in% levels(enthesis$name2))) {
    message(paste("Running Slingshot using reduction:", TRAJECTORY_REDUCTION))
    enthesis <- SCP::RunSlingshot(srt = enthesis, group.by = "name2",
                                  start_node = TRAJECTORY_START_NODE,
                                  end_node = TRAJECTORY_END_NODES,
                                  reduction = TRAJECTORY_REDUCTION,
                                  show_plot = TRUE) 
    
    lineage_tool_name <- names(enthesis@tools)[grep("^Slingshot$", names(enthesis@tools))][[1]]
    actual_lineage_names <- names(slot(enthesis@tools[[lineage_tool_name]], "lineages"))
    
    if (length(actual_lineage_names) > 0) {
      # Select first 2 distinct lineages for plotting if more exist (original had Lineage1, Lineage3)
      lineages_to_plot <- actual_lineage_names[1:min(2, length(actual_lineage_names))]
      lineage_plot_colors <- RColorBrewer::brewer.pal(max(3, length(lineages_to_plot)), "Set1")[1:length(lineages_to_plot)]
      
      pdf(ENTHESIS_LINEAGE_PLOT_PDF, width = 8, height = 7)
      print(scplotter::CellDimPlot(enthesis, group_by = "name2", reduction = TRAJECTORY_REDUCTION,
                                   lineages = lineages_to_plot,
                                   lineages_trim = c(0.01, 0.99), lineages_whiskers = FALSE, lineages_span = 0.1,
                                   palcolor = name2_cols, pt_size =1.0,
                                   lineages_palcolor = lineage_plot_colors,
                                   theme_light = FALSE))
      dev.off()
      
      # Dynamic Gene Features
      DefaultAssay(enthesis) <- 'RNA' # Ensure RNA assay for expression
      if (length(actual_lineage_names) > 0) { # Run on all identified lineages
        enthesis <- SCP::RunDynamicFeatures(srt = enthesis, lineages = actual_lineage_names,
                                            n_candidates = DYN_FEAT_N_CANDIDATES)
        
        dynamic_genes_to_plot_scp <- c('Tnmd',"Wif1",'Scrg1','Tnn') 
        dynamic_genes_present_scp <- dynamic_genes_to_plot_scp[dynamic_genes_to_plot_scp %in% rownames(enthesis)]
        
        if (length(dynamic_genes_present_scp) > 0 && length(lineages_to_plot) > 0) {
          pdf(ENTHESIS_DYNFEAT_PLOT_PDF, width = 5, height = 3 * length(dynamic_genes_present_scp))
          print(SCP::DynamicPlot(srt = enthesis, lineages = lineages_to_plot, group.by = "name2",
                                 line_palcolor = lineage_plot_colors, point_palcolor = name2_cols,
                                 features = dynamic_genes_present_scp, ncol = 1,
                                 compare_lineages = TRUE, compare_features = FALSE))
          dev.off()
        }
        # Dynamic Heatmap 
        if (length(actual_lineage_names) >= 2) {
          ht_plot_obj <- SCP::DynamicHeatmap(
            srt = enthesis, lineages = actual_lineage_names[1:2],
            use_fitted = TRUE, GO_simplify = FALSE, # db_update = TRUE (check SCP docs for this param)
            split_method = "mfuzz", n_split = 4, # reverse_ht = actual_lineage_names[1],
            species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE,
            heatmap_palette = "RdBu", cell_annotation = "name2", cell_annotation_palcolor = name2_cols,
            height = 4, width = 1.5
          )
          pdf(ENTHESIS_DYNHEATMAP_GO_PDF, width = 15, height = 10)
          print(ht_plot_obj$plot) # Assuming $plot contains the ggplot object or ComplexHeatmap object
          dev.off()
        }
      }
    } else { message("Slingshot did not identify any lineages.") }
  } else { message("Start or end nodes for trajectory not found in 'name2' levels. Skipping Slingshot.") }
} else { message(paste("Skipping Slingshot: 'name2' column or reduction '", TRAJECTORY_REDUCTION, "' not found.", sep=""))}

# --- 10.3 GO Enrichment with SCP::RunEnrichment ---
message("--- Part 10.3: GO Enrichment with SCP ---")
if ("name2" %in% colnames(enthesis@meta.data)) {
  # DefaultAssay(enthesis) <- "SCT" # Or 'integrated'
  # enthesis <- SCP::RunDEtest(srt = enthesis, group_by = "name2", assay = DefaultAssay(enthesis),
  #                            fc.threshold = log2(1.5), p_val_adj.threshold = 0.05, only.pos = FALSE) # Example thresholds
  
  if ("DEtest_name2" %in% names(enthesis@tools)) { 
    enthesis <- SCP::RunEnrichment(
      srt = enthesis, group_by = "name2", db = "GO_BP", species = "Mus_musculus",
      DE_threshold_gene = "avg_log2FC > log2(1.5) & p_val_adj < 0.05", 
      enrich_method = "hypergeometric" # Default is hypergeometric
    )
    
    pdf(file.path(RESULTS_DIR_DEV, "scp_name2_go_enrichment_comparison.pdf"), width = 10, height = 8)
    print(SCP::EnrichmentPlot(srt = enthesis, group_by = "name2", plot_type = "comparison", topn_terms = 5))
    dev.off()
  } else {
    message("DE test results for 'name2' not found. Skipping SCP::RunEnrichment.")
  }
}

# --- 10.4 Subset 'Chondrogenesis' Population and Analyze ---
message("--- Part 10.4: Subsetting and Analyzing 'Chondrogenesis' Population ---")
if ("name2" %in% colnames(enthesis@meta.data)) {
  chondro_cell_types <- c('Enthesis chondrocytes','Enthesis progenitors')
  if (all(chondro_cell_types %in% levels(enthesis$name2))) {
    Idents(enthesis) <- 'name2'
    chondrogenesis <- subset(enthesis, idents = chondro_cell_types)
    chondrogenesis@project.name <- "Chondrogenesis_Subset"
    
    DefaultAssay(chondrogenesis) <- 'RNA'
    chondrogenesis <- NormalizeData(chondrogenesis, verbose = FALSE)
    chondrogenesis <- ScaleData(chondrogenesis, features = rownames(chondrogenesis),
                                vars.to.regress = c('percent.mt','nCount_RNA','nFeature_RNA'), verbose = FALSE)
    saveRDS(chondrogenesis, CHONDROGENESIS_SUBSET_RDS)
    # chondrogenesis <- readRDS(CHONDROGENESIS_SUBSET_RDS)
    
    # FeatureStatPlot for specific genes by time
    pdf(file.path(RESULTS_DIR_DEV, "chondro_subset_featurestat_time.pdf"), width=8, height=4)
    print(scplotter::FeatureStatPlot(chondrogenesis, assay = 'RNA', features = c("Scx", "Sox9",'Tnn'),
                                     group_by = "time_point", plot_type = "violin", add_bg = TRUE, add_box = TRUE, stack = TRUE))
    dev.off()
    
  } else {
    message("Chondrogenesis cell types not found in 'name2'. Skipping subsetting.")
  }
}

# --- 10.5 Trajectory Analysis with GeneTrajectory Package ---
message("--- Part 10.5: Trajectory Analysis with GeneTrajectory ---")
if (exists("chondrogenesis") && inherits(chondrogenesis, "Seurat")) {
  tryCatch({
    library(GeneTrajectory)
    library(plot3D) 
    library(SeuratWrappers) # For RunALRA
    
    DefaultAssay(chondrogenesis) <- "RNA"
    chondrogenesis <- FindVariableFeatures(chondrogenesis, nfeatures = GENETRAJ_VAR_FEATURES)
    
    # Filter variable features
    var_genes_gt <- chondrogenesis@assays[["RNA"]]@var.features
    patterns_to_exclude <- c('^Rp[sl]', "\\b\\w+Rik\\b", '^Hist', '^Gm', '^Ccl', '^Cxcl', '^Cdc',
                             '^Hbb','^Hba','^Ccn','^Hmg','^Ifi','^C1q','^C4','^Pi1','^Il',
                             '^Atp','^Ndu','^Tmem','^Tnf','^Cdk','H2')
    var_genes_gt_filtered <- var_genes_gt[!sapply(var_genes_gt, function(g) any(sapply(patterns_to_exclude, function(p) grepl(p, g, ignore.case=TRUE))))]
    
    expr_percent_gt <- Matrix::rowSums(chondrogenesis[["RNA"]]@data[var_genes_gt_filtered, ] > 0) / ncol(chondrogenesis)
    genes_for_gt <- var_genes_gt_filtered[expr_percent_gt > GENETRAJ_EXPR_PERCENT_MIN & expr_percent_gt < GENETRAJ_EXPR_PERCENT_MAX]
    
    if (length(genes_for_gt) > 10) { # Need enough genes
      chondrogenesis <- GeneTrajectory::RunDM(chondrogenesis, K = 20) # Uses "RNA" assay, HVGs
      
      # Setup Python env for GeneTrajectory
      # reticulate::py_config() # Check config
      use_python(GENETRAJ_PYTHON_PATH, required = TRUE)
      # use_virtualenv(GENETRAJ_VENV_PATH, required=TRUE) # If it's a venv
      
      cal_ot_mat_from_numpy_fn <- reticulate::import('gene_trajectory.compute_gene_distance_cmd')$cal_ot_mat_from_numpy
      
      cell_graph_dist_gt <- GetGraphDistance(chondrogenesis, K = GENETRAJ_K, dims = 1:10) # Uses DM reduction
      cg_output_gt <- CoarseGrain(chondrogenesis, cell_graph_dist_gt, genes_for_gt, N = GENETRAJ_META_CELLS_N)
      
      gene_dist_mat_gt <- cal_ot_mat_from_numpy_fn(ot_cost = cg_output_gt[["graph.dist"]],
                                                   gene_expr = cg_output_gt[["gene.expression"]],
                                                   num_iter_max = 50000L, show_progress_bar = TRUE)
      rownames(gene_dist_mat_gt) <- colnames(gene_dist_mat_gt) <- cg_output_gt[["features"]]
      
      gene_embedding_gt <- GetGeneEmbedding(gene_dist_mat_gt, K = GENETRAJ_K)$diffu.emb
      gene_trajectory_gt <- ExtractGeneTrajectory(gene_embedding_gt, gene_dist_mat_gt, N = GENETRAJ_TRAJECTORIES_N,
                                                  t.list = c(1,1,2), K = GENETRAJ_K) # t.list might need tuning
      
      # Visualization
      pdf(GENETRAJ_SCATTER3D_PDF, width = 7, height = 7)
      par(mar = c(1.5,1.5,1.5,1.5))
      plot3D::scatter3D(gene_embedding_gt[,3], gene_embedding_gt[,1], gene_embedding_gt[,2],
                        bty = "b2", colvar = as.integer(as.factor(gene_trajectory_gt$selected)),
                        main = "Gene Trajectory (GeneTrajectory)", pch = 19, cex = 1.2, theta = 60, phi = -10,
                        col = RColorBrewer::brewer.pal(GENETRAJ_TRAJECTORIES_N, "Set1"))
      dev.off()
      
      # ALRA and GeneBinScore
      chondrogenesis <- RunALRA(chondrogenesis) # Imputes data, stores in "alra" assay
      chondrogenesis <- AddGeneBinScore(chondrogenesis, gene_trajectory_gt,
                                        N.bin = 5, trajectories = 1, assay = "alra", reverse = c(TRUE,TRUE))
      
      pdf(GENETRAJ_FEATUREPLOT_PDF, width = 15, height = 3)
      print(FeaturePlot(chondrogenesis, pt.size = 0.05, features = paste0("Trajectory",1,"_genes", 1:5),
                        reduction = 'umap', ncol = 5, order = TRUE) &
              scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n=9, name="RdYlBu"))) & NoLegend() & NoAxes())
      dev.off()
    } else { message("Not enough variable genes meeting expression criteria for GeneTrajectory.") }
  }, error = function(e) { message(paste("GeneTrajectory analysis failed:", e$message))})
} else { message("Skipping GeneTrajectory: 'chondrogenesis' subset not available.") }

# --- 10.6 Data Export ---
message("--- Part 10.6: Data Export ---")
if (exists("enthesis") && inherits(enthesis, "Seurat")) {
  # Ensure time factor levels are set as desired for export
  enthesis$time_point <- factor(x=enthesis$time_point, levels = intersect(c('e15','p01','p07','p14','p28'), unique(enthesis$time_point))) # p01 in original, maps to p07?
  
  # Save as H5Seurat
  SaveH5Seurat(enthesis, filename = ENTHESIS_H5SEURAT_OUT, overwrite = TRUE)
  # Convert to AnnData (RNA assay, data slot)
  Convert(ENTHESIS_H5SEURAT_OUT, dest = "h5ad", assay = 'RNA', main_layer="data", overwrite = TRUE,
          outfile = ENTHESIS_H5AD_OUT_RNA) # Specify output file for Convert
  
  # Ensure Python env for sceasy is active if needed
  # use_python(...)
  sceasy::convertFormat(enthesis, from="seurat", to="anndata", assay="RNA",
                        main_layer="counts", drop_single_values=FALSE, outFile = ENTHESIS_H5AD_OUT_SCVI) # Original path
  message(paste("Exported 'enthesis' to H5Seurat and AnnData formats in:", RESULTS_DIR_DEV, "and", dirname(ENTHESIS_H5AD_OUT_SCVI)))
}

# --- 10.7 CytoTRACE2 Analysis (on 'tnn' object) ---
message("--- Part 10.7: CytoTRACE2 Analysis ---")
# This section seems to operate on a different object 'tnn'.
# Ensure 'tnn' is loaded and is a Seurat object.
if (exists("tnn") && inherits(tnn, "Seurat")) {
  tryCatch({
    library(CytoTRACE2)
    library(ggpubr) # For stat_compare_means
    
    cytotrace_result_tnn <- cytotrace2(tnn, assay_ epigenomics = "RNA", # is_seurat = TRUE is default
                                       slot_type = "counts", species = 'mouse', seed = 1234)
    
    # Add results to tnn object if desired
    tnn$CytoTRACE2_Score_obj <- cytotrace_result_tnn$CytoTRACE2_Score
    tnn$CytoTRACE2_Potency_obj <- cytotrace_result_tnn$CytoTRACE2_Potency
    
    # Plotting (original used plotData, which returns a list of plots)
    annotation_tnn <- data.frame(phenotype = tnn$time_point) # Assuming tnn has 'time_point'
    rownames(annotation_tnn) <- colnames(tnn)
    
    plots_ct_tnn <- plotData(cytotrace2_result = cytotrace_result_tnn, is_seurat = TRUE, annotation = annotation_tnn)
    
    # Custom boxplot comparison
    df_ct_tnn <- plots_ct_tnn$CytoTRACE2_Boxplot_byPheno$data
    # Define comparisons (example, adapt to your 'time_point' levels in 'tnn')
    # my_comparisons_tnn <- list(c("e15", "p07"), c("p07", "p14")) # Example
    my_comparisons_tnn <- list() 
    if(length(unique(df_ct_tnn$Phenotype)) >= 2) {
      my_comparisons_tnn <- combn(unique(as.character(df_ct_tnn$Phenotype)), 2, simplify = FALSE)
    }
    
    
    pdf(CYTOTRACE_PLOTS_PDF, height = 5, width = 7)
    if(length(my_comparisons_tnn) > 0) {
      print(ggplot(df_ct_tnn, aes(Phenotype, CytoTRACE2_Score, fill = Phenotype)) +
              geom_boxplot(alpha = 0.8, outlier.shape = NA) +
              stat_compare_means(comparisons = my_comparisons_tnn, method = 't.test', label="p.signif", size=4) +
              scale_fill_manual(values = sample_cols) + 
              theme_bw() + labs(title="CytoTRACE Score by Time (tnn object)", y = "CytoTRACE Score") +
              theme(axis.text.x = element_text(angle=45, hjust=1)))
    } else {
      print(ggplot(df_ct_tnn, aes(Phenotype, CytoTRACE2_Score, fill = Phenotype)) +
              geom_boxplot(alpha = 0.8, outlier.shape = NA) +
              scale_fill_manual(values = sample_cols) + 
              theme_bw() + labs(title="CytoTRACE Score by Time (tnn object)", y = "CytoTRACE Score") +
              theme(axis.text.x = element_text(angle=45, hjust=1)))
    }
    dev.off()
    
  }, error = function(e) { message(paste("CytoTRACE2 analysis on 'tnn' failed:", e$message))})
} else { message("Skipping CytoTRACE2 analysis: 'tnn' object not found or not a Seurat object.") }


# --- 10.8 Final Cleanup ---
message("--- Part 10.8: Final Cleanup ---")
gc()

sessionInfo()
message(paste0("Advanced developmental analysis script execution finished at ", Sys.time(), "! Results are in: ", RESULTS_DIR_DEV))