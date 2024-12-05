library(here)
library(SeuratDisk)
library(BiocManager)
library(Seurat)
library(SeuratObject)
library(dplyr)


#DEG ANALSIS
setwd("C:/Users/mahsa/Documents/data/interim/rds_VIS")


# Check if running interactively
if (interactive()) {
  # If in RStudio or interactive mode, define default arguments
  args <- c("ClusteredVisium_LogNorm.rds")  # Specify a default input file here
} else {
  # Reads the argument provided at the command line
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) == 0) {
  stop("An argument must be supplied (input file).")
}

input_file <- args[1]
cat("Using input file:", input_file, "\n")

data <- readRDS(here("data", "interim", "rds_VIS",
                     
                     args[1]))



#If running Differential Expression on the clusters of Visium data, the
# default assay needs to be set to 'Spatial' (instead of 'RNA')
DefaultAssay(data) <- "Spatial"
class(data)
# Setting group of interest and idents
group <- "disease"
Idents(data) <- "seurat_clusters"

# Change variable 'top' to desired number of differentially expressed genes
top <- 500
data <- LogNormalize(data)

marker <- list()

for (cluster in levels(data$seurat_clusters)) {
  # Subset data to current cluster
  cluster_subset <- subset(data, subset = seurat_clusters == cluster)
  
  # Get list of disease groups (visium -> HC, NL, LS; scRNAseq -> HC, AD)
  group_id <- levels(as.factor(cluster_subset$disease))
  
  # Set ident to "disease"
  Idents(cluster_subset) <- "disease"
  
  # Checks if cluster is present in only one disease condition
  if (length(group_id) < 2) {
    # Skip differential expression if cluster contains less than 3 cells
    if (length(WhichCells(cluster_subset,
                          idents = as.character(group_id[1]))) < 3 &&
        length(WhichCells(cluster_subset,
                          idents = as.character(group_id[2]))) < 3) {
      next
    }
    
    # Find markers for cluster with >= 3 cells or spots (but only present in
    # one disease)
    x <- FindMarkers(data, slot = "data", ident.1 = cluster,
                     only.pos = TRUE)
    
    # Select top 500 differentially expressed genes by adjusted p-value
    x <- x %>% slice_min(p_val_adj, n = top)
    
    # Add cluster id as column
    x$cluster <- cluster
    
    # Append data frame to list 'marker'
    marker[[length(marker) + 1]] <- x
    rm(x)
    next
  } else if (length(WhichCells(cluster_subset,
                               idents = as.character(group_id[1]))) < 3 ||
             length(WhichCells(cluster_subset,
                               idents = as.character(group_id[2]))) < 3) {
    
    # Find markers for cluster with >= 3 cells or spots in one disease only
    x <- FindMarkers(data, slot = "data", ident.1 = cluster,
                     only.pos = TRUE)
    
    # Select top 500 differentially expressed genes by adjusted p-value
    x <- x %>% slice_min(p_val_adj, n = top)
    
    # Add cluster id as column
    x$cluster <- cluster
    
    # Append data frame to list 'marker'
    marker[[length(marker) + 1]] <- x
    rm(x)
    next
  }
  
  # Set assay to default assay, since FindConservedMarkers() will overwrite the
  # default assay with "RNA" by default (only necessary for Visium DE testing).
  assay <- DefaultAssay(data)
  
  # Finds conserved markers across groups for the currently selected cluster.
  # Note: FindConservedMarkers only works if there are at least 3 cells/spots
  # in each disease group, hence the previous if-else shenanigans and the use of
  # FindMarkers if this condition is not met
  x <- FindConservedMarkers(data, assay = assay, slot = "data",
                            ident.1 = cluster, grouping.var = group,
                            only.pos = TRUE)
  
  # Keep the top 500 markers (by combined p-value = 'minimump_p_pvalue')
  x <- x %>% slice_min(minimump_p_val, n = top)
  
  # Adds the cluster id as column to the DE output data frame, that can later be
  # used to view the top distinguishing genes for this cluster
  x$cluster <- cluster
  
  # Adds the data frame to the marker list
  marker[[length(marker) + 1]] <- x
  rm(x)
}


saveRDS(marker,
        file = here("DEGs_IntegratedRO_LogNorm_1stLvl.rds"),
        compress = TRUE)


############################################################################################################
#assign cell type

scdata.integrated <- readRDS(
  here("ClusteredVisium_LogNorm.rds"))

# Restores the saved list of data.frames (DE output).
Ref.markers <- readRDS(here( "DEGs_IntegratedRO_LogNorm_1stLvl.rds"))

DefaultAssay(cluster0) <- "integrated"

while (TRUE) {
  cluster <- readline(prompt = "Enter cluster number or x to exit the loop: ")
  if (cluster == "x") {
    break
  } else {
    df <- Ref.markers[[as.numeric(cluster)+1]]
    df$DEG <- row.names(df)
    tib <- as_tibble(df)
    print(tib, n = nrow(tib))
  }
}


# Renames cluster identities to their putative cell types.
scdata.integrated$cluster_name <- recode(scdata.integrated$seurat_clusters,
                                         "0" = "KC", "1" = "Lymphoid",
                                         "2" = "Myeloid", "3" = "FB",
                                         "4" = "MEL", "5" = "VEC",
                                         "6" = "SMC", "7" = "MC",
                                         "8" = "LEC", "9" = "PDC")





saveRDS(
  scdata.integrated,
  file = here("final_clustreing.rds"),
  compress = TRUE)

metadata_df <-scdata.integrated@meta.data
write.csv(metadata_df, file = "fianal_meta_data.csv", row.names = TRUE)
