# ADDITIONAL INFORMATION PROVIDED HERE

# Seurat v5.1.0.
# SeuratObject v5.0.2.


# Load Libraries ----
# Close libaries
if("Seurat" %in% (.packages())){
  detach("package:Seurat", unload=TRUE) 
}

if("SeuratObject" %in% (.packages())){
  detach("package:SeuratObject", unload=TRUE) 
}

library(clustree)
library(cowplot)
library(dittoSeq)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(glmGamPoi)
library(harmony)
library(hdf5r)
library(matrixStats)
library(metafolio)
library(patchwork)
# library(presto)
library(RColorBrewer)
library(readxl)
library(sctransform)
# library(SeuratDisk)
library(Seurat)
# library(SeuratData)
library(writexl)


print("!!! LIBRARIES LOADED SUCCESSFULLY. !!!")

# STEP 1: Set working directory & sample information. ----
setwd("E:/MSc_EMC_project/Main_project/02_BatchCorr_Niches_outs/filtereddata/")
data_path = "E:/MSc_EMC_project/Main_project/SR_outs/"

# SAMPLE INFORMATION
info_samples.xlsx <- read_excel("E:/MSc_EMC_project/Main_project/Main File_snRNAseq cases_v3xlsx.xlsx",
                                sheet = "Sheet1", range = "B1:E139")
colnames(info_samples.xlsx) <- c("PA_number", "Diagnosis", "ZIS", "Gender")

print("!!! STEP 1 FINISHED. !!!")

# STEP 2: Load data & find HVGs. ----

# Load data
filtered_combined_samples <- readRDS("E:/MSc_EMC_project/Main_Project/01_QC_outs/Visium_samples_filtered_combined_QC_5removed.RDS")
# Normalize 
filtered_combined_samples <- NormalizeData(filtered_combined_samples, normalization.method = 'LogNormalize',
                                  scale.factor = 10000, verbose = F)

# Find highly variable genes (HVGs)
# Create empty HVGs list
hvgs_list <- NULL

# Set identities to original identities.
Idents(filtered_combined_samples) <- "orig.ident"

# # Set image keys to corresponding sample
# for(x in 1:length(names(filtered_combined_samples@images))){
#   filtered_combined_samples@images[[x]]@key <- names(filtered_combined_samples@images)[[x]]
# }

# Find HGVs for each sample
for(sample in unique(filtered_combined_samples$orig.ident)) {
  # For each sample, subset it with the corresponding sample identity.
  x <- subset(filtered_combined_samples, idents = sample)
  
  # Find variable features for each sample.
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 3000)
  
  # Add to list under sample name.
  hvgs_list[[sample]] <- x@assays$Spatial@meta.data$var.features
}

# Unlist and turn into character.
hvgs_list <- unlist(hvgs_list)

# Create HVGs table, sort by decreasing gene frequency.
hvgs_list <- table(hvgs_list) %>%
  sort(decreasing = T)

# Get names of top 3000 genes.
gene_selection <- hvgs_list[1:3000] %>% 
  names()

print("!!! STEP 2 FINISHED. !!!")

# STEP 3: Dimensionality reduction & visualization. ----
# PCA prep.
filtered_combined_samples <- ScaleData(filtered_combined_samples, verbose = F)
# Perform PCA using top 3000 highly variable genes.
filtered_combined_samples <- RunPCA(filtered_combined_samples, verbose = F,
                           npcs = 30, features = gene_selection)

plot_1 <- ElbowPlot(filtered_combined_samples, ndims = 20) + ggtitle("PCA elbow plot")
ggsave("01_ElbowPlot_PCA.jpeg",
       width = 5)

print("!!! STEP 3 FINISHED. !!!")

# STEP 4: Harmony integration & clustering. ----

## 4.1. Perform Harmony integration & clustering. ----
filtered_combined_samples <- RunHarmony(filtered_combined_samples,
                               group.by.vars = "orig.ident", 
                               plot_convergence = T,
                               assay.use = "Spatial", max_iter = 20)

# Elbow plot post-Harmony integration.
plot_1 <- ElbowPlot(filtered_combined_samples, ndims = 20, reduction = "harmony") +
  ggtitle("Harmony integration elbow plot")
elbow_harmony <- plot_1
ggsave("02_ElbowPlot_Harmony.jpeg", plot_1,
       width = 5)

# Simple clustering based on Harmony dimensions
filtered_combined_samples <- RunUMAP(filtered_combined_samples, reduction = "harmony",
                            dims = 1:15)
filtered_combined_samples <- FindNeighbors(filtered_combined_samples, 
                                  reduction = "harmony", dims = 1:15)
filtered_combined_samples <- identity(filtered_combined_samples)

# Dimplots
plot_1 <- DimPlot(filtered_combined_samples, group.by = "diag", label = T, pt.size = 1,
                  reduction = "umap")
ggsave("03_DimPlot_umap_diag.jpeg", plot_1,
       height = 10, width = 12)

plot_1 <- DimPlot(filtered_combined_samples, label = T, pt.size = 1,
                  reduction = "umap")
ggsave("04_DimPlot_umap.jpeg", plot_1,
       height = 10, width = 12)

## 4.2. Clustering ----
# Set resolutions
resolutions <- seq(0.1, 1.0, 0.1)
# Find clusters in resolutions range
filtered_combined_samples <- FindClusters(filtered_combined_samples, resolution = resolutions,
                                 verbose = F)

# Optimize clustering
# Obtain cell distances
cell_dists <- dist(filtered_combined_samples@reductions$harmony@cell.embeddings,
                   method = "euclidean")

cluster_info <- filtered_combined_samples@meta.data[,grepl(paste0(DefaultAssay(filtered_combined_samples),"_snn_res"),
                                                  colnames(filtered_combined_samples@meta.data))]
cluster_info <- mutate_all(cluster_info, as.character)
cluster_info <- mutate_all(cluster_info, as.numeric)

# Apply function over cluster_info columns (0.1-1), score for each clustering resolution.
silhouette_res <- apply(cluster_info, 2, function(x) {
  # Silhouette width for each resolution.
  si <- cluster::silhouette(x, cell_dists)
  # If not NA, calculate mean.
  if (!any(is.na(si[, 'sil_width']))) {
    mean(si[, 'sil_width'])
    # Else, return NA.
  } else {
    NA
  }
})

# Turn into data frame.
silhouette_res <- as.data.frame(silhouette_res)
# Replace column names.
silhouette_res$res <- gsub("Spatial_snn_res.","", rownames(silhouette_res))
silhouette_res$res <- as.numeric(silhouette_res$res)

# Plot results
plot_1 <- ggplot(silhouette_res, aes(res,silhouette_res))+
  geom_point()+
  theme_classic()
ggsave("05_silhouette_ress.jpeg", plot_1,
       width = 5 , height = 5)

# Remove items redundant in future analysis for memory space.
rm(cell_dists)
rm(x)
rm(cluster_info)
rm(hvgs_list)

print("!!! STEP 4 FINISHED !!!")



# STEP 5: Loop through resolutions & plot. ----
## 5.1. Find anchors ----
# # Load reference data obtained from KPMP atlas.
# ref_data <-  LoadH5Seurat("E:/MSc_EMC_project/Main_Project/12be772d-49bd-4abc-b08e-76684c2ef1f2_KPMPAtlas_PREMIER_062921.h5Seurat")
# # Update the Seurat object
# ref_data <- UpdateSeuratObject(ref_data)
# 
# # Save as RDS to obtain data more quickly in the future.
# saveRDS(ref_data, "refdata_anchors.RDS")

# Load the reference data.
ref_data <- readRDS("E:/MSc_EMC_project/Main_project/02_BatchCorr_Niches_outs/refdata_anchors.RDS")
# Normalize
ref_data <- NormalizeData(ref_data,
                          normalization.method = 'LogNormalize', 
                          scale.factor = 10000, verbose = F)
DefaultAssay(ref_data) <- "RNA"

# Source: https://satijalab.org/seurat/archive/v4.3/spatial_vignette
# Predict the major celltype for each bead.
anchors <- FindTransferAnchors(reference = ref_data, query = filtered_combined_samples,
                               normalization.method = "LogNormalize", npcs = 50)

# Standard Seurat label transfer workflow
predictions <- TransferData(anchorset = anchors, refdata = ref_data$celltype, prediction.assay = TRUE,
                            weight.reduction = filtered_combined_samples[["pca"]], dims = 1:30)

filtered_combined_samples[["predictions"]] <- predictions
DefaultAssay(filtered_combined_samples) <- "predictions"

prediction_results <- t(as.data.frame(predictions@data))

# Row names correspond to the cell names of the Seurat object 
# Column names correspond to the metadata variables
filtered_combined_samples <- AddMetaData(filtered_combined_samples, metadata = prediction_results)

### Plots for transfer anchors.

for (anchor_no in 1:length(colnames(prediction_results))){
  print(paste0("Current transfer anchor: ", colnames(prediction_results)[anchor_no]))
  plot_list <- list()
  
  for (image in names(filtered_combined_samples@images)){
    plot <- SpatialFeaturePlot(filtered_combined_samples, features = colnames(prediction_results)[anchor_no], images = image,
                               alpha = c(0.6, 1), crop = F) + ggtitle(paste0(image))
    plot_list[[image]] <- plot
    
    ggsave(plot = plot, filename = paste0("transfer_anchors/", colnames(prediction_results)[anchor_no], "/", image, ".jpeg"),
           width = 10, height = 10)
  }
  ggsave(wrap_plots(plot_list, ncol = 5), filename = paste0("transfer_anchors/", colnames(prediction_results)[anchor_no], "_wrapped.jpeg"),
         width = 13, height = 7)
  plot_list <- list()
}


## 5.2. Loop through resolutions ----
DefaultAssay(filtered_combined_samples) <- "Spatial"

# Create empty plot list and set x to 1.
plot_list <- list()
x <- 1
# Try different resolutions: 0.3-0.7
for (res in seq(0.3, 0.7, 0.1)) {
  print(paste0("Current resolution: ", res))
  
  # Find clusters at the set resolution
  filtered_combined_samples <- FindClusters(filtered_combined_samples, resolution = res, verbose = F)
  
  plot_list[[x]] <- DimPlot(filtered_combined_samples, label = T, label.size = 8) +
    NoLegend() + labs(title = paste0("UMAP at resolution ", res))

  
  # Find marker genes
  all_markers <- FindAllMarkers(filtered_combined_samples, min.pct = 0.3, logfc.threshold = 0.3,
                                assay = "Spatial", only.pos = T)
  # Obtain genes with positive log2FC & with adj p-values below 0.01
  all_markers <- all_markers[all_markers$p_val_adj < 0.01, ]
  
  # Obtain top 10 marker genes per cluster.
  top10 <- all_markers %>% 
    # Group marker genes by cluster number.
    group_by(cluster) %>% 
    # Obtain top 10 genes per cluster based on avg log2FC in cluster.
    top_n(n = 10, wt = avg_log2FC) %>% 
    # Extract gene names.
    pull(gene)
  
  filtered_combined_samples <- ScaleData(filtered_combined_samples, features = unique(all_markers$gene))
  seurat_cluster_colors <- gg_color_hue(length(levels(filtered_combined_samples$seurat_clusters)))
  
  # Plot heatmaps
  heatmap <- DoHeatmap(filtered_combined_samples,
                       features = top10, group.by = "seurat_clusters",
                       slot = "scale.data", label = TRUE) +
    theme(axis.text.y = element_text(size = 7))
  
  # Save grouped heatmap.
  ggsave(filename = paste0("niche_res/", res, "_Grouped_CM_subcluster_allmarkers_heatmap.jpeg"),
         plot = heatmap, width = 20, height = 20, limitsize = F)
  
  x <- x + 1
  
}

resolutions_plot_list <- plot_list

plot_1 <- wrap_plots(plot_list)
ggsave(filename = "niche_res/CM_subcluster_umap_seurat_cluster_resolutions.jpeg", plot_1,
       width = 20 , height = 20)

# Visualize clusters using clustree.
x <- clustree(filtered_combined_samples) + labs(title = "Clustree all samples")
clustree_plot <- clustree(filtered_combined_samples) + labs(title = "    Clustree all samples")

ggsave("06_clustree_filtered_combined_samples.jpeg", x,
       width = 10, height = 10)

# Save after Harmony integration and clustering.
saveRDS(filtered_combined_samples, "Visium_filtered_combined_samples_harmony_corrected.RDS")

print("!!! STEP 5 FINISHED !!!")





# STEP 6: Use chosen resolution & plot. ----
filtered_combined_samples <- readRDS("./Visium_filtered_combined_samples_harmony_corrected.RDS")

# Find clusters using chosen resolution.
filtered_combined_samples <- FindClusters(filtered_combined_samples, resolution = 0.4)
# Find marker genes
all_markers <- FindAllMarkers(filtered_combined_samples, min.pct = 0.3, logfc.threshold = 0.3,
                              assay = "Spatial")
# Save as Excel file.
write_xlsx(all_markers, path = "filtered_combined_samples_AllMarkers.xlsx")

# Obtain genes with positive log2FC & with adj p-values below 0.01
all_markers <- all_markers[all_markers$avg_log2FC > 0 & all_markers$p_val_adj < 0.01, ]
# Save as Excel file.
write_xlsx(all_markers, path = "filtered_combined_samples_PosMarkers_p-val_below_0.01.xlsx")


## 6.1. Dotplot/Heatmap TOP 3 GENES per cluster. ----
# Obtain top 3 genes per cluster based on avg log2FC in cluster.
gene_clusters <- group_by(all_markers, cluster)
top3_gene_clusters <- top_n(gene_clusters, n = 3, wt = avg_log2FC)


# Dotplots of top3 genes per cluster.
plot <- DotPlot(filtered_combined_samples, features = top3_gene_clusters$gene,
        group.by = 'seurat_clusters') + labs(title = paste0("Top 3 genes per cluster")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
top3_dotplot <- plot

# Save as Excel file & save dotplot.
write_xlsx(top3_gene_clusters, path = "dotplots/top3_Dotplot.xlsx")
ggsave(filename = paste0("dotplots/top3_Dotplot.jpeg"),
       plot = plot, limitsize = F, width = 12)

# Heatmap of top3 genes per cluster.
heatmap <- DoHeatmap(filtered_combined_samples,
                     features = top3_gene_clusters$gene, group.by = "seurat_clusters", size = 4,
                     slot = "scale.data", label = TRUE) + 
  labs(title = paste0("Top 3 genes per cluster")) 
ggsave(filename = paste0("dotplots/top3_heatmap.jpeg"),
       plot = heatmap, limitsize = F, width = 12)

## 6.2. Dotplot/Heatmap TOP 10 GENES per cluster. ----
# Obtain top 10 genes per cluster based on avg log2FC in cluster.
top10_gene_clusters <- top_n(gene_clusters, n = 10, wt = avg_log2FC)

# Dotplots of top10 genes per cluster.
plot <- DotPlot(filtered_combined_samples, features = top10_gene_clusters$gene,
                group.by = 'seurat_clusters') + labs(title = paste0("Top 10 genes per cluster")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
top10_dotplot <- plot

# Save as Excel file & save dotplot.
write_xlsx(top10_gene_clusters, path = "dotplots/top10_Dotplot.xlsx")
ggsave(filename = paste0("dotplots/top10_Dotplot.jpeg"),
       plot = plot, limitsize = F, width = 25)

# Heatmap of top3 genes per cluster.
heatmap <- DoHeatmap(filtered_combined_samples,
                     features = top10_gene_clusters$gene, group.by = "seurat_clusters", size = 4,
                     slot = "scale.data", label = TRUE) + 
  labs(title = paste0("Top 10 genes per cluster")) 
ggsave(filename = paste0("dotplots/top10_heatmap.jpeg"),
       plot = heatmap, limitsize = F, width = 12, height = 15)


print("!!! STEP 6 FINISHED !!!")



## Dotplots questionable niches ----

plot <- DotPlot(filtered_combined_samples, features = c('COL1A1', 'POSTN', 'COL3A1', 'FN1'),
                group.by = 'seurat_clusters') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_1.jpeg"),
       plot = plot, limitsize = F, width = 10)

plot <- VlnPlot(filtered_combined_samples, features = c('COL1A1', 'POSTN', 'COL3A1', 'FN1'),
                group.by = 'seurat_clusters') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_1_VLN.jpeg"),
       plot = plot, limitsize = F, width = 10, height = 10)

plot <- DotPlot(filtered_combined_samples, features = c('COL1A1', 'POSTN', 'COL3A1', 'FN1'),
                group.by = 'diag') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_1_diag.jpeg"),
       plot = plot, limitsize = F, width = 10)

plot <- VlnPlot(filtered_combined_samples, features = c('COL1A1', 'POSTN', 'COL3A1', 'FN1'),
                group.by = 'diag') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_1_VLN_diag.jpeg"),
       plot = plot, limitsize = F, width = 10, height = 10)

#//////////////////////////////////////////////////////////////////////////////////////////////
plot <- DotPlot(filtered_combined_samples, features = c('SPP1', 'CD14', 'CD16', 'ARG1',
                                                        'HIF1A', 'ITGAM', 'CD96', 'F13A1',
                                                        'MRC1', 'CD163', 'MSR1', 'ITGAX',
                                                        'TPRM2', 'S100A8', 'S100A9'),
                group.by = 'seurat_clusters') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_fibrotic-immune.jpeg"),
       plot = plot, limitsize = F, width = 10)

plot <- VlnPlot(filtered_combined_samples, features = c('SPP1', 'CD14', 'CD16', 'ARG1',
                                                        'HIF1A', 'ITGAM', 'CD96', 'F13A1',
                                                        'MRC1', 'CD163', 'MSR1', 'ITGAX',
                                                        'TPRM2', 'S100A8', 'S100A9'),
                group.by = 'seurat_clusters') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_fibrotic-immune_VLN.jpeg"),
       plot = plot, limitsize = F, width = 15, height = 15)

plot <- DotPlot(filtered_combined_samples, features = c('SPP1', 'CD14', 'CD16', 'ARG1',
                                                        'HIF1A', 'ITGAM', 'CD96', 'F13A1',
                                                        'MRC1', 'CD163', 'MSR1', 'ITGAX',
                                                        'TPRM2', 'S100A8', 'S100A9'),
                group.by = 'diag') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_fibrotic-immune_diag.jpeg"),
       plot = plot, limitsize = F, width = 10)

plot <- VlnPlot(filtered_combined_samples, features = c('SPP1', 'CD14', 'CD16', 'ARG1',
                                                        'HIF1A', 'ITGAM', 'CD96', 'F13A1',
                                                        'MRC1', 'CD163', 'MSR1', 'ITGAX',
                                                        'TPRM2', 'S100A8', 'S100A9'),
                group.by = 'diag') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_fibrotic-immune_VLN_diag.jpeg"),
       plot = plot, limitsize = F, width = 15, height = 15)


#//////////////////////////////////////////////////////////////////////////////////////////////
plot <- DotPlot(filtered_combined_samples, features = c('CD34', 'PECAM1', 'PTPRB', 'EMCN', 'FLT1'),
                group.by = 'seurat_clusters') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_endothelial.jpeg"),
       plot = plot, limitsize = F, width = 10)

plot <- VlnPlot(filtered_combined_samples, features = c('CD34', 'PECAM1', 'PTPRB', 'EMCN', 'FLT1'),
                group.by = 'seurat_clusters') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_endothelial_VLN.jpeg"),
       plot = plot, limitsize = F, width = 10, height = 10)

plot <- DotPlot(filtered_combined_samples, features = c('CD34', 'PECAM1', 'PTPRB', 'EMCN', 'FLT1'),
                group.by = 'diag') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_endothelial_diag.jpeg"),
       plot = plot, limitsize = F, width = 10)

plot <- VlnPlot(filtered_combined_samples, features = c('CD34', 'PECAM1', 'PTPRB', 'EMCN', 'FLT1'),
                group.by = 'diag') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)
ggsave(filename = paste0("dotplots/niches_endothelial_VLN_diag.jpeg"),
       plot = plot, limitsize = F, width = 10, height = 10)


# STEP 7: Barplots. ----

# Barplot per DIAGNOSIS.
barplot_1 <- dittoBarPlot(object = filtered_combined_samples,
                          var = "seurat_clusters",
                          group.by = "diag",
                          scale = 'count',
                          main = "Composition of clusters per diagnosis (count)",
                          legend.title = "Cluster no.") 

barplot_2 <- dittoBarPlot(object = filtered_combined_samples,
                          var = "seurat_clusters",
                          group.by = "diag",
                          scale = 'percent',
                          main = "Composition of clusters per diagnosis (%)",
                          legend.title = "Cluster no.") 

ggsave("07_Barplots_bydiag.jpeg" , wrap_plots(barplot_1, barplot_2),
       width = 20, height = 10)


# Barplot per SAMPLE.
barplot_3 <- dittoBarPlot(object = filtered_combined_samples,
                          var = "seurat_clusters",
                          group.by = "orig.ident",
                          scale = 'count',
                          main = "Composition of clusters per sample (count)",
                          legend.title = "Cluster no.")

barplot_4 <- dittoBarPlot(object = filtered_combined_samples,
                          var = "seurat_clusters",
                          group.by = "orig.ident",
                          scale = 'percent',
                          main = "Composition of clusters per sample (%)",
                          legend.title = "Cluster no.") 

ggsave("08_Barplots_bysample.jpeg" , wrap_plots(barplot_3, barplot_4),
       width = 20, height = 10)

# Save all barplots.
ggsave("09_Barplots.jpeg", wrap_plots(barplot_1, barplot_2,barplot_3, barplot_4),
       height = 20, width = 20)

ggsave("09.5_Barplots.jpeg", wrap_plots(barplot_2,barplot_4),
       height = 15, width = 20)

print("!!! STEP 7 FINISHED !!!")


# STEP 8: SpatialDimPlots. ----

p <- distinct(filtered_combined_samples@meta.data, orig.ident, diag)

plot_list <- list()
x <- 1

for(sample in names(filtered_combined_samples@images)){
  
  d <- p$diag[which(p$orig.ident == sample)]
  
  plot_list[[x]] <- SpatialDimPlot(filtered_combined_samples, group.by = "seurat_clusters", label = TRUE, repel = T,
                                   images = sample, crop = F, pt.size.factor =  1.6, label.size = 4) + 
    NoLegend() + labs(title = paste0(sample, " (", d, ")"))
  
  x <- x + 1
}

# Wrap and save Dimplot per Visium sample.
ggsave("10.1_SpatialDimPlots.jpeg", wrap_plots(plot_list[1:2]),
       width = 15, height = 9)
ggsave("10.2_SpatialDimPlots.jpeg", wrap_plots(plot_list[3:5]),
       width = 25, height = 9)
ggsave("10.3_SpatialDimPlots.jpeg", wrap_plots(plot_list[6:7]),
       width = 15, height = 9)
ggsave("10.4_SpatialDimPlots.jpeg", wrap_plots(plot_list[8:10]),
       width = 25, height = 9)
# Save all.
ggsave("10.5_SpatialDimPlots_ALL.jpeg", wrap_plots(plot_list, ncol = 3), 
       height = 30, width = 20)

print("!!! STEP 8 FINISHED !!!")


# STEP 9: Volcano plots ----
## 9.1. TCMR vs IF/TA ----
Idents(filtered_combined_samples) <- "diag"
TCMRvsIFTA <- FindMarkers(filtered_combined_samples,
                          ident.1 = c('aTCMR1B','aTCMR2B'), 
                          ident.2 = "IFTA")
TCMRvsIFTA <- TCMRvsIFTA[TCMRvsIFTA$p_val_adj <= 0.05,]
TCMRvsIFTA$p_val_adj[TCMRvsIFTA$p_val_adj == 0] <- 1e-300

# Add genes as separate column.
TCMRvsIFTA$gene <- rownames(TCMRvsIFTA)

TCMRvsIFTA$up_down <- "NA"
TCMRvsIFTA$up_down[TCMRvsIFTA$avg_log2FC > 0.6] <- "UP"
TCMRvsIFTA$up_down[TCMRvsIFTA$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(TCMRvsIFTA, (avg_log2FC > 3 | avg_log2FC < -3  | p_val_adj < (1/(10^290)) & p_val_adj != 1e-300)) 

volcano_TCMRvsIFTA <- ggplot(TCMRvsIFTA, aes(x = avg_log2FC, y = -log10(p_val_adj),
                       colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 5, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "aTCMR vs IF/TA", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf, colour = "black") +
  guides(colour = guide_legend(override.aes = aes(label = ""))) 

# Save as Excel file & save volcano plot.
write_xlsx(TCMRvsIFTA, path = "volcano_plots/01_TCMRvsIFTA.xlsx")
ggsave("volcano_plots/01_TCMRvsIFTA.jpeg", plot = volcano_TCMRvsIFTA,
       height = 10, width = 10)
 

## 9.2. AMR vs IF/TA ----
AMRvsIFTA <- FindMarkers(filtered_combined_samples,
                         ident.1 = 'aAMR, C4d+',
                         ident.2 = 'IFTA')
AMRvsIFTA <- AMRvsIFTA[AMRvsIFTA$p_val_adj <= 0.05,]
AMRvsIFTA$p_val_adj[AMRvsIFTA$p_val_adj == 0] <- 1e-300
# Add genes as separate column.
AMRvsIFTA$gene <- rownames(AMRvsIFTA)

AMRvsIFTA$up_down <- "NA"
AMRvsIFTA$up_down[AMRvsIFTA$avg_log2FC > 0.6] <- "UP"
AMRvsIFTA$up_down[AMRvsIFTA$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(AMRvsIFTA, (avg_log2FC > 6 | avg_log2FC < -5.5 | p_val_adj < (1/(10^200)) | 
                                 avg_log2FC > 3 & p_val_adj < (1/(10^150)) | avg_log2FC > 2.8))

volcano_AMRvsIFTA <- ggplot(AMRvsIFTA, aes(x = avg_log2FC, y = -log10(p_val_adj),
                      colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 5, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "aAMR vs IF/TA", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf, colour = "black") +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(AMRvsIFTA, path = "volcano_plots/02_AMRvsIFTA.xlsx")
ggsave("volcano_plots/02_AMRvsIFTA.jpeg", plot = volcano_AMRvsIFTA,
       height = 10, width = 10)

## 9.3. TCMR vs AMR ----
TCMRvsAMR <- FindMarkers(filtered_combined_samples,
                         ident.1 = c('aTCMR1B', 'aTCMR2B'), ident.2 = 'aAMR, C4d+')
TCMRvsAMR <- TCMRvsAMR[TCMRvsAMR$p_val_adj <= 0.05,]
TCMRvsAMR$p_val_adj[TCMRvsAMR$p_val_adj == 0] <- 1e-300
# Add genes as separate column.
TCMRvsAMR$gene <- rownames(TCMRvsAMR)

TCMRvsAMR$up_down <- "NA"
TCMRvsAMR$up_down[TCMRvsAMR$avg_log2FC > 0.6] <- "UP"
TCMRvsAMR$up_down[TCMRvsAMR$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(TCMRvsAMR, (avg_log2FC > 7 | avg_log2FC < -3 | p_val_adj < (1/(10^235))))

volcano_TCMRvsAMR <- ggplot(TCMRvsAMR, aes(x = avg_log2FC, y = -log10(p_val_adj),
                      colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 12, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "aTCMR vs aAMR", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf, colour = "black") +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(TCMRvsAMR, path = "volcano_plots/03_TCMRvsAMR.xlsx")
ggsave("volcano_plots/03_TCMRvsAMR.jpeg", plot = volcano_TCMRvsAMR,
       height = 10, width = 10)



# Thesis figures ----
## Figure xx - Volcano plots ----
# Fig_VolcanoPlots <- ggarrange(ggarrange((volcano_TCMRvsIFTA + NoLegend()), (volcano_AMRvsIFTA + NoLegend()), ncol = 2, labels = c("A", "B")),
#                               ggarrange(volcano_TCMRvsAMR, labels = "C"), 
#                               font.label = list(size = 20))
# # OPTION 1
# Fig_VolcanoPlots <- ggarrange((volcano_IFTAvsTCMR + NoLegend()), (volcano_IFTAvsAMR + NoLegend()), volcano_TCMRvsAMR,
#                                         labels = c("A", "B", "C"), 
#                               font.label = list(size = 20), ncol = 3)
# ggsave(plot = Fig_VolcanoPlots, filename = "Thesis_Fig_2.1.jpeg",
#        height = 8, width = 25)

# OPTION 2
Fig_VolcanoPlots2 <- ggarrange(ggarrange((volcano_TCMRvsIFTA + NoLegend()), (volcano_AMRvsIFTA + NoLegend()), labels = c("A", "B")),
                              ggarrange(volcano_TCMRvsAMR, labels = "C"), nrow = 2, 
                              font.label = list(size = 20))
ggsave(plot = Fig_VolcanoPlots2, filename = "Thesis_Fig_2.223_2.jpeg",
       height = 11, width = 13)

## Figure xx - volcano plot tables ----
library(gt)

# TCMR vs IF/TA
TCMRvsIFTA_table_up <- top_n(TCMRvsIFTA, n = 5, wt = avg_log2FC)
TCMRvsIFTA_table_up <- TCMRvsIFTA_table_up[order(TCMRvsIFTA_table_up$avg_log2FC, decreasing = TRUE), ]

TCMRvsIFTA_table_down <- top_n(TCMRvsIFTA, n = 5, wt = -avg_log2FC)
TCMRvsIFTA_table_down <- TCMRvsIFTA_table_down[order(TCMRvsIFTA_table_down$avg_log2FC), ]

TCMRvsIFTA_table <- rbind(TCMRvsIFTA_table_up, TCMRvsIFTA_table_down)
TCMRvsIFTA_table <- TCMRvsIFTA_table[,2:6]
TCMRvsIFTA_table <- TCMRvsIFTA_table[,-2:-3]
colnames(TCMRvsIFTA_table) <- c("Log2FC", "Adjusted p-value", "Gene")

# Create tables
gt(TCMRvsIFTA_table) |> tab_header(title = md("**Upregulated and downregulated genes in TCMR vs IF/TA**"),
                                   subtitle = "Top 5 highest and lowest Log2FC are presented") |>
  tab_row_group(label = md("**Upregulated**"), rows = 1:5) |>
  tab_row_group(label = md("**Downregulated**"), rows = 6:10) |>
  cols_move_to_start(columns = "Gene") |>
  tab_options(table.font.names = "Cambria", table.font.size = 12) |>
  fmt_number(decimals = , columns = "Log2FC") |>
  tab_options(data_row.padding = pct(1)) |>
  gt::gtsave(filename = "Table3_volcano_output.png")


# AMR vs IF/TA
AMRvsIFTA_table_up <- top_n(AMRvsIFTA, n = 5, wt = avg_log2FC)
AMRvsIFTA_table_up <- AMRvsIFTA_table_up[order(AMRvsIFTA_table_up$avg_log2FC, decreasing = TRUE), ]

AMRvsIFTA_table_down <- top_n(AMRvsIFTA, n = 5, wt = -avg_log2FC)
AMRvsIFTA_table_down <- AMRvsIFTA_table_down[order(AMRvsIFTA_table_down$avg_log2FC), ]

AMRvsIFTA_table <- rbind(AMRvsIFTA_table_up, AMRvsIFTA_table_down)
AMRvsIFTA_table <- AMRvsIFTA_table[,2:6]
AMRvsIFTA_table <- AMRvsIFTA_table[,-2:-3]
colnames(AMRvsIFTA_table) <- c("Log2FC", "Adjusted p-value", "Gene")

gt(AMRvsIFTA_table) |> tab_header(title = md("**Upregulated and downregulated genes in AMR vs IF/TA**"),
                                   subtitle = "Top 5 highest and lowest Log2FC are presented") |>
  tab_row_group(label = md("**Upregulated**"), rows = 1:5) |>
  tab_row_group(label = md("**Downregulated**"), rows = 6:10) |>
  cols_move_to_start(columns = "Gene") |>
  tab_options(table.font.names = "Cambria", table.font.size = 12) |>
  fmt_number(decimals = , columns = "Log2FC") |>
  tab_options(data_row.padding = pct(1)) |>
  gt::gtsave(filename = "Table4_volcano_output.png")


# TCMR vs AMR
TCMRvsAMR_table_up <- top_n(TCMRvsAMR, n = 5, wt = avg_log2FC)
TCMRvsAMR_table_up <- TCMRvsAMR_table_up[order(TCMRvsAMR_table_up$avg_log2FC, decreasing = TRUE), ]

TCMRvsAMR_table_down <- top_n(TCMRvsAMR, n = 5, wt = -avg_log2FC)
TCMRvsAMR_table_down <- TCMRvsAMR_table_down[order(TCMRvsAMR_table_down$avg_log2FC), ]

TCMRvsAMR_table <- rbind(TCMRvsAMR_table_up, TCMRvsAMR_table_down)
TCMRvsAMR_table <- TCMRvsAMR_table[,2:6]
TCMRvsAMR_table <- TCMRvsAMR_table[,-2:-3]
colnames(TCMRvsAMR_table) <- c("Log2FC", "Adjusted p-value", "Gene")

gt(TCMRvsAMR_table) |> tab_header(title = md("**Upregulated and downregulated genes in TCMR vs AMR**"),
                                  subtitle = "Top 5 highest and lowest Log2FC are presented") |>
  tab_row_group(label = md("**Upregulated**"), rows = 1:5) |>
  tab_row_group(label = md("**Downregulated**"), rows = 6:10) |>
  cols_move_to_start(columns = "Gene") |>
  tab_options(table.font.names = "Cambria", table.font.size = 12) |>
  fmt_number(decimals = , columns = "Log2FC") |>
  tab_options(data_row.padding = pct(1)) |>
  gt::gtsave(filename = "Table5_volcano_output.png")



## Figure xx - Barplot ----
barplot_figure <- ggarrange(barplot_2, barplot_4, labels = c("A", "B"), font.label = list(size = 20))
ggsave(plot = barplot_figure, filename = "Thesis_Fig_barplot.jpeg",
       height = 6, width = 10)

## Figure xx - Harmony, UMAP, dotplot, clustree ----
# Variables
# elbow_harmony
# resolutions_plot_list
# top3_dotplot
# top10_dotplot
# clustree_plot

clustering_figure <- ggarrange(elbow_harmony,
                               ggarrange(resolutions_plot_list[[2]], top3_dotplot, nrow = 2),
                               clustree_plot, ncol = 3)
# OPTION 1
clustering_figure <- ggarrange(resolutions_plot_list[[2]],
                               ggarrange(elbow_harmony, top3_dotplot, nrow = 2, labels = c("B", "C")),
                               clustree_plot, ncol = 3, labels = c("A", "", "D"),
                               font.label = list(size = 20))
ggsave(plot = clustering_figure, filename = "Thesis_Fig_clustering.1.jpeg",
       height = 10 , width = 25)

# OPTION 2
clustering_figure <- ggarrange(elbow_harmony,
                               ggarrange(resolutions_plot_list[[2]], top3_dotplot, nrow = 2, labels = c("B", "C")),
                               clustree_plot, ncol = 3, labels = c("A", "", "D"),
                               font.label = list(size = 20))
ggsave(plot = clustering_figure, filename = "Thesis_Fig_clustering.2.jpeg",
       height = 10 , width = 25)

# OPTION 3
clustering_figure <- ggarrange(
  ggarrange(elbow_harmony,resolutions_plot_list[[2]], nrow = 2, labels = c("A", "B")),
  ggarrange(, top3_dotplot, nrow = 2, labels = c("C", "D")),
  ggarrange(barplot_2, barplot_4, nrow = 2, labels = c("E", "F")),
  ncol = 3,
  font.label = list(size = 20))
ggsave(plot = clustering_figure, filename = "Thesis_Fig_clustering.3.jpeg",
       height = 15 , width = 25)

# OPTION 4
clustering_figure <- ggarrange(
  ggarrange(elbow_harmony,clustree_plot, nrow = 2, labels = c("A", "B"), font.label = list(size = 20)),
  resolutions_plot_list[[2]], labels = c("", "C"),
  ggarrange(barplot_2, barplot_4, nrow = 2, labels = c("D", "E"), font.label = list(size = 20)),
  ncol = 3,
  font.label = list(size = 20))
ggsave(plot = clustering_figure, filename = "Thesis_Fig_clustering.4.jpeg",
       height = 18 , width = 25)

# OPTION 5
# clustering_figure <- ggarrange(
#   ggarrange(elbow_harmony,clustree_plot, nrow = 2, labels = c("A", "B"), font.label = list(size = 15)),
#   resolutions_plot_list[[2]], labels = c("", "C"),
#   ggarrange(barplot_2, barplot_4, nrow = 2, labels = c("D", "E"), font.label = list(size = 15)),
#   ncol = 3,
#   font.label = list(size = 15))
# ggsave(plot = clustering_figure, filename = "Thesis_Fig_clustering.5.jpeg",
#        height = 18 , width = 25)



## Figure xx - dotplot ----
library(gt)

cluster_no <- c(as.character(0:10))
annotation <- c("Tubular: loop of Henle",
                "Tubular: Proximal tubules (1)",
                "Tubular: Distal tubules",
                "Fibrotic",
                "Tubular: Proximal tubules (2)",
                "Mixed",
                "Collecting duct",
                "Immune niche",
                "Glomeruli",
                "Vascular niche",
                "Renal capsule")

annotation_table <- data.frame(cluster_no, annotation)
colnames(annotation_table) <- c("Cluster", "Annotation")

gt(annotation_table) |> 
  tab_options(table.font.names = "Cambria", table.font.size = 25) |>
  gtsave("test.png")

dotplot_fig <- ggarrange(
  ggarrange(top3_dotplot, ggdraw() + draw_image("test.png"), labels = c("A", "C"), widths = c(3,1),
            font.label = list(size = 20)),
  top10_dotplot, nrow = 2, labels = c("", "B"),
  font.label = list(size = 20))

ggsave(plot = dotplot_fig, filename = "Thesis_Fig3_dotplots_annotation.jpeg",
       width = 20, height = 10)



# Removed from analysis / OLD ----
## 9.4. IFTA vs [AMR and HLAi] ----
IFTAvsAMR_HLAi <- FindMarkers(filtered_combined_samples,
                              ident.1 = "IFTA", ident.2 = c('aAMR, C4d+','HLAi'))
IFTAvsAMR_HLAi <- IFTAvsAMR_HLAi[IFTAvsAMR_HLAi$p_val_adj != 0 & IFTAvsAMR_HLAi$p_val_adj <= 0.05,]
# Add genes as separate column.
IFTAvsAMR_HLAi$gene <- rownames(IFTAvsAMR_HLAi)

IFTAvsAMR_HLAi$up_down <- "NA"
IFTAvsAMR_HLAi$up_down[IFTAvsAMR_HLAi$avg_log2FC > 0.6] <- "UP"
IFTAvsAMR_HLAi$up_down[IFTAvsAMR_HLAi$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(IFTAvsAMR_HLAi, (avg_log2FC > 4 | avg_log2FC < -3 | p_val_adj < (1/(10^200))))

to_save <- ggplot(IFTAvsAMR_HLAi, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                      colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 6, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "IFTA vs [AMR and HLAi]", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf) +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(IFTAvsAMR_HLAi, path = "volcano_plots/02_IFTAvsAMR_HLAi.xlsx")
ggsave("volcano_plots/04_IFTAvsAMR_HLAi.jpeg", plot = to_save,
       width = 10)

## 9.5. IFTA vs HLAi ----
IFTAvsHLAi <- FindMarkers(filtered_combined_samples,
                          ident.1 = 'IFTA', ident.2 = 'HLAi')
IFTAvsHLAi <- IFTAvsHLAi[IFTAvsHLAi$p_val_adj != 0 & IFTAvsHLAi$p_val_adj <= 0.05,]
# Add genes as separate column.
IFTAvsHLAi$gene <- rownames(IFTAvsHLAi)

IFTAvsHLAi$up_down <- "NA"
IFTAvsHLAi$up_down[IFTAvsHLAi$avg_log2FC > 0.6] <- "UP"
IFTAvsHLAi$up_down[IFTAvsHLAi$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(IFTAvsHLAi, (avg_log2FC > 3.5 | avg_log2FC < -3 | p_val_adj < (1/(10^150))))

to_save <- ggplot(IFTAvsHLAi, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                  colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 6, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "IFTA vs HLAi", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf) +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(IFTAvsHLAi, path = "volcano_plots/04_IFTAvsHLAi.xlsx")
ggsave("volcano_plots/05_IFTAvsHLAi.jpeg", plot = to_save,
       width = 10)

## 9.6. TCMR vs [AMR and HLAi] ----
TCMRvsAMR_HLAi <- FindMarkers(filtered_combined_samples,
                              ident.1 = c("aTCMR1B", "aTCMR2B"), ident.2 = c('aAMR, C4d+','HLAi'))
TCMRvsAMR_HLAi <- TCMRvsAMR_HLAi[TCMRvsAMR_HLAi$p_val_adj != 0 & TCMRvsAMR_HLAi$p_val_adj <= 0.05,]
# Add genes as separate column.
TCMRvsAMR_HLAi$gene <- rownames(TCMRvsAMR_HLAi)

TCMRvsAMR_HLAi$up_down <- "NA"
TCMRvsAMR_HLAi$up_down[TCMRvsAMR_HLAi$avg_log2FC > 0.6] <- "UP"
TCMRvsAMR_HLAi$up_down[TCMRvsAMR_HLAi$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(TCMRvsAMR_HLAi, (avg_log2FC > 5 | avg_log2FC < -3 | p_val_adj < (1/(10^235))))

to_save <- ggplot(TCMRvsAMR_HLAi, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                      colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 6, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "TCMR vs [AMR and HLAi]", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf) +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(TCMRvsAMR_HLAi, path = "volcano_plots/05_TCMRvsAMR_HLAi.xlsx")
ggsave("volcano_plots/06_TCMRvsAMR_HLAi.jpeg", plot = to_save,
       width = 10)

## 9.7. TCMR vs HLAi ----
TCMRvsHLAi <- FindMarkers(filtered_combined_samples,
                          ident.1 = c('aTCMR1B', 'aTCMR2B'), ident.2 = 'HLAi')
TCMRvsHLAi <- TCMRvsHLAi[TCMRvsHLAi$p_val_adj != 0 & TCMRvsHLAi$p_val_adj <= 0.05,]

# Add genes as separate column.
TCMRvsHLAi$gene <- rownames(TCMRvsHLAi)

TCMRvsHLAi$up_down <- "NA"
TCMRvsHLAi$up_down[TCMRvsHLAi$avg_log2FC > 0.6] <- "UP"
TCMRvsHLAi$up_down[TCMRvsHLAi$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(TCMRvsHLAi, (avg_log2FC > 5 | avg_log2FC < -2 | p_val_adj < (1/(10^200))))

to_save <- ggplot(TCMRvsHLAi, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                  colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 6, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "TCMR vs HLAi", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf) +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(TCMRvsHLAi, path = "volcano_plots/07_TCMRvsHLAi.xlsx")
ggsave("volcano_plots/07_TCMRvsHLAi.jpeg", plot = to_save,
       width = 10)

print("!!! STEP 9 FINISHED. !!!")



## STEP 9: Volcano plots ----
## 9.1. IFTA vs TCMR ----
Idents(filtered_combined_samples) <- "diag"
IFTAvsTCMR <- FindMarkers(filtered_combined_samples,
                          ident.1 = "IFTA", ident.2 = c('aTCMR1B','aTCMR2B'))
IFTAvsTCMR <- IFTAvsTCMR[IFTAvsTCMR$p_val_adj <= 0.05,]
# Add genes as separate column.
IFTAvsTCMR$gene <- rownames(IFTAvsTCMR)

IFTAvsTCMR$up_down <- "NA"
IFTAvsTCMR$up_down[IFTAvsTCMR$avg_log2FC > 0.6] <- "UP"
IFTAvsTCMR$up_down[IFTAvsTCMR$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(IFTAvsTCMR, (avg_log2FC > 3 | avg_log2FC < -3 | p_val_adj < (1/(10^290))))

volcano_IFTAvsTCMR <- ggplot(IFTAvsTCMR, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                             colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 5, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "IFTA vs TCMR", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf, colour = "black") +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(IFTAvsTCMR, path = "volcano_plots/01_IFTAvsTCMR.xlsx")
ggsave("volcano_plots/01_IFTAvsTCMR.jpeg", plot = volcano_IFTAvsTCMR,
       height = 10, width = 10)


## 9.2. IFTA vs AMR ----
IFTAvsAMR <- FindMarkers(filtered_combined_samples,
                         ident.1 = "IFTA",
                         ident.2 = 'aAMR, C4d+')
IFTAvsAMR <- IFTAvsAMR[IFTAvsAMR$p_val_adj <= 0.05,]
# Add genes as separate column.
IFTAvsAMR$gene <- rownames(IFTAvsAMR)

IFTAvsAMR$up_down <- "NA"
IFTAvsAMR$up_down[IFTAvsAMR$avg_log2FC > 0.6] <- "UP"
IFTAvsAMR$up_down[IFTAvsAMR$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(IFTAvsAMR, (avg_log2FC > 6 | avg_log2FC < -2 | p_val_adj < (1/(10^200)) | 
                                 avg_log2FC > 3 & p_val_adj < (1/(10^150))))

volcano_IFTAvsAMR <- ggplot(IFTAvsAMR, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                           colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 8, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "IFTA vs AMR", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf, colour = "black") +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(IFTAvsAMR, path = "volcano_plots/02_IFTAvsAMR.xlsx")
ggsave("volcano_plots/02_IFTAvsAMR.jpeg", plot = volcano_IFTAvsAMR,
       height = 10, width = 10)

## 9.3. TCMR vs AMR ----
TCMRvsAMR <- FindMarkers(filtered_combined_samples,
                         ident.1 = c('aTCMR1B', 'aTCMR2B'), ident.2 = 'aAMR, C4d+')
TCMRvsAMR <- TCMRvsAMR[TCMRvsAMR$p_val_adj <= 0.05,]
# Add genes as separate column.
TCMRvsAMR$gene <- rownames(TCMRvsAMR)

TCMRvsAMR$up_down <- "NA"
TCMRvsAMR$up_down[TCMRvsAMR$avg_log2FC > 0.6] <- "UP"
TCMRvsAMR$up_down[TCMRvsAMR$avg_log2FC < -0.6] <- "DOWN"

# Label genes of interest.
interest <- filter(TCMRvsAMR, (avg_log2FC > 7 | avg_log2FC < -3 | p_val_adj < (1/(10^235))))

volcano_TCMRvsAMR <- ggplot(TCMRvsAMR, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                           colour = up_down, label = gene)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') +
  annotate("text", x = 12, y = 8, label = 'p = 0.05', col = "gray47") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  scale_color_manual(values = c("cornflowerblue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(title = "TCMR vs AMR", color = "") +
  geom_text_repel(data = interest, max.overlaps = Inf, colour = "black") +
  guides(colour = guide_legend(override.aes = aes(label = "")))

# Save as Excel file & save volcano plot.
write_xlsx(TCMRvsAMR, path = "volcano_plots/03_TCMRvsAMR.xlsx")
ggsave("volcano_plots/06_TCMRvsAMR.jpeg", plot = volcano_TCMRvsAMR,
       height = 10, width = 10)


