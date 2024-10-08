# Load libaries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(clustermole)
library(writexl)

# Set directory & load data ----
setwd("E:/MSc_EMC_project/Main_project/02_BatchCorr_Niches_outs/filtereddata/")
filtered_combined_samples <- readRDS("./Visium_filtered_combined_samples_harmony_corrected.RDS")


# Find clusters ----
filtered_combined_samples <- FindClusters(filtered_combined_samples, resolution = 0.4)
# Find marker genes
all_markers <- FindAllMarkers(filtered_combined_samples, min.pct = 0.3, logfc.threshold = 0.3,
                              assay = "Spatial")
# Obtain genes with positive log2FC & with adj p-values below 0.01
all_markers <- all_markers[all_markers$avg_log2FC > 0 & all_markers$p_val_adj < 0.01, ]
# Group per cluster and obtain top 10 marker genes per cluster.
gene_clusters <- group_by(all_markers, cluster)
top10_gene_clusters <- top_n(gene_clusters, n = 10, wt = avg_log2FC)

# Obtain annotated marker genes ----
markers <- clustermole_markers(species = "hs")
kidney_markers <- filter(markers, organ == 'Kidney' & species == 'Human')
# Set gene names to uppercase to match top10_gene_clusters
kidney_markers$gene <- toupper(kidney_markers$gene)

# Obtain present markers ----
annotated_genes <- inner_join(top10_gene_clusters, kidney_markers, by = "gene")
results <- bind_cols(annotated_genes$cluster, annotated_genes$celltype,
                     annotated_genes$gene_original, annotated_genes$db)
colnames(results) <- c("Cluster", "Cell type", "Gene", "Database")

write_xlsx(results, path = "Cluster_annotation.xlsx")
