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


library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(glmGamPoi)
library(gprofiler2)
library(readxl)
library(sctransform)
library(Seurat)
library(stringr)
require(svglite)
library(tidyverse)
# library(Vennerable)
library(writexl)

print("!!! LIBRARIES LOADED SUCCESSFULLY. !!!")

# STEP 1: Set working directory & sample information. ----
setwd("E:/MSc_EMC_project/Main_project/03_GO_analysis_outs/")
starting_directory <- getwd()

# SAMPLE INFORMATION
info_samples.xlsx <- read_excel("E:/MSc_EMC_project/Main_project/Main File_snRNAseq cases_v3xlsx.xlsx",
                                sheet = "Sheet1", range = "B1:E139")
colnames(info_samples.xlsx) <- c("PA_number", "Diagnosis", "ZIS", "Gender")

print("!!! STEP 1 FINISHED. !!!")

# STEP 2: Set parameters. ----

# Define diagnoses.
IFTA <- "IFTA"
TCMR <- c("aTCMR1B", "aTCMR2B")
AMR <- "aAMR, C4d+"

print("!!! STEP 2 FINISHED. !!!")

# STEP 3: Load data & cluster ----
filtered_combined_samples <- readRDS("E:/MSc_EMC_project/Main_project/02_BatchCorr_Niches_outs/filtereddata/Visium_filtered_combined_samples_harmony_corrected.RDS")
# Find clusters using chosen resolution.
filtered_combined_samples <- FindClusters(filtered_combined_samples, resolution = 0.4)

# Add column with cluster no. and diagnosis.
filtered_combined_samples$clusterno_diag <- paste(filtered_combined_samples$seurat_clusters,
                                                filtered_combined_samples$diag, sep = "_")
Idents(filtered_combined_samples) <- "clusterno_diag"

print("!!! STEP 3 FINISHED. !!!")

# STEP 4: GO Analysis. ----
## 4.1. TCMR vs IF/TA ----
setwd(starting_directory)
# dir.create("TCMRvsIFTA")
setwd(paste0(starting_directory, "/TCMRvsIFTA"))

Idents(filtered_combined_samples) <- 'diag'
TCMRvsIFTA <- FindMarkers(filtered_combined_samples,
                          test.use = "MAST",
                          min.pct = 0.25, logfc.threshold = 0.3,
                          ident.1 = TCMR, 
                          ident.2 = IFTA)

# Sort by avg_log2FC.
TCMRvsIFTA <- arrange(TCMRvsIFTA, -avg_log2FC)
# Add gene names as separate column.
TCMRvsIFTA$gene <- rownames(TCMRvsIFTA)
# Save as Excel file & RDS.
write_xlsx(TCMRvsIFTA, path = paste0("DE_FindMarkers_gProfiler_TCMRvsIFTA.xlsx"))
saveRDS(TCMRvsIFTA, paste0("DE_FindMarkers_gProfiler.RDS"))

### UPREGULATED ----
# Obtain upregulated gene names.
upregulated <- TCMRvsIFTA$gene[TCMRvsIFTA$avg_log2FC > 0 & 
                                 TCMRvsIFTA$p_val_adj < 0.01]
# GO analysis using upregulated genes.
GO_upregulated <- gost(query = upregulated, 
                       organism = "hsapiens", ordered_query = T,
                       correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_upregulated <- GO_upregulated$result
# Save results.
write_xlsx(results_GO_upregulated, path = paste0("TCMRvsIFTA_results_GO_upregulated.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
TCMRvsIFTA_GO_up <- BP_GO_up_filtered
# Obtain top 3/10 GOs by lowest p_value.
GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
TCMRvsIFTA_GO_top3  <- top_n(BP_GO_up_filtered , n = -5, wt = p_value)
# Obtain top 3/10 term names.
top10_term_names_GO_up <- unique(GO_up_top10$term_name)
TCMRvsIFTA_GO_top3_terms  <- unique(TCMRvsIFTA_GO_top3$term_name)


# Plot results
TCMRvsIFTA_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                              y = factor(source), size = precision, fill = -log10(p_value),
                                              color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                             subtitle = "TCMR vs IF/TA") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = TCMRvsIFTA_up_plot, filename = "TCMRvsIFTA_GO_up.jpeg")


### DOWNREGULATED ----
# Obtain downregulated gene names.
downregulated <- TCMRvsIFTA$gene[TCMRvsIFTA$avg_log2FC < 0 & 
                                   TCMRvsIFTA$p_val_adj < 0.01]
# GO analysis using downregulated genes.
GO_downregulated <- gost(query = downregulated, 
                         organism = "hsapiens", ordered_query = T,
                         correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_downregulated <- GO_downregulated$result
# Save results.
write_xlsx(results_GO_downregulated, path = paste0("TCMRvsIFTA_results_GO_downregulated.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_down_filtered <- filter(results_GO_downregulated, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_down_filtered <- dplyr::filter(GO_down_filtered, source == "GO:BP")
# Obtain top 10 GOs by lowest p_value.
GO_down_top10 <- top_n(BP_GO_down_filtered, n = -10, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_down <- unique(GO_down_top10$term_name)

# Plot results
TCMRvsIFTA_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                  y = factor(source), size = precision, fill = -log10(p_value),
                                                  color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                             subtitle = "TCMR vs IF/TA") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = TCMRvsIFTA_down_plot, filename = "TCMRvsIFTA_GO_down.jpeg")

### ALL ----
# Obtain all gene names.
all_genes <- TCMRvsIFTA$gene[TCMRvsIFTA$p_val_adj < 0.01]
# GO analysis using all genes.
GO_allgenes <- gost(query = all_genes, 
                    organism = "hsapiens", ordered_query = T,
                    correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_all <- GO_allgenes$result
# Save results.
write_xlsx(results_GO_all, path = paste0("TCMRvsIFTA_results_GO_all.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_all_filtered <- filter(results_GO_all, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_all_filtered <- dplyr::filter(GO_all_filtered, source == "GO:BP")
# Obtain top 10 GOs by lowest p_value.
GO_all_top10 <- top_n(BP_GO_all_filtered, n = -10, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_all <- unique(GO_all_top10$term_name)

# Plot results
TCMRvsIFTA_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                 y = factor(source), size = precision, fill = -log10(p_value),
                                                 color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                             subtitle = "TCMR vs IF/TA") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = TCMRvsIFTA_down_plot, filename = "TCMRvsIFTA_GO_all.jpeg")

### Save results ----
TCMRvsIFTA_GO_results <- NULL
TCMRvsIFTA_GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
saveRDS(TCMRvsIFTA_GO_results,
        file = paste0("gProfiler.RDS"))

# Plots top 10 of upregulated - downregulated - all
# UP
GO_upregulated$result <- BP_GO_up_filtered
plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                     subtitle = "TCMR vs IF/TA")
publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_upregulated.jpeg"))

# DOWN
GO_downregulated$result <- BP_GO_down_filtered
plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                       subtitle = "TCMR vs IF/TA")
publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_downregulated.jpeg"))
#ALL
GO_allgenes$result <- BP_GO_all_filtered
plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                  subtitle = "TCMR vs IF/TA")
publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_all.jpeg"))

## 4.2. AMR vs IF/TA ----
setwd(starting_directory)
# dir.create("AMRvsIFTA")
setwd(paste0(starting_directory, "/AMRvsIFTA"))

Idents(filtered_combined_samples) <- 'diag'
AMRvsIFTA <- FindMarkers(filtered_combined_samples,
                         test.use = "MAST",
                         min.pct = 0.25, logfc.threshold = 0.3,
                         ident.1 = AMR, 
                         ident.2 = IFTA)

# Sort by avg_log2FC.
AMRvsIFTA <- arrange(AMRvsIFTA, -avg_log2FC)
# Add genes as separate column.
AMRvsIFTA$gene <- rownames(AMRvsIFTA)
# Save as Excel file & RDS.
write_xlsx(AMRvsIFTA, path = paste0("DE_FindMarkers_gProfiler_AMRvsIFTA.xlsx"))
saveRDS(AMRvsIFTA, paste0("DE_FindMarkers_gProfiler.RDS"))

### UPREGULATED ----
# Obtain upregulated gene names.
upregulated <- AMRvsIFTA$gene[AMRvsIFTA$avg_log2FC > 0 & 
                                AMRvsIFTA$p_val_adj < 0.01]
# GO analysis using upregulated genes.
GO_upregulated <- gost(query = upregulated, 
                       organism = "hsapiens", ordered_query = T,
                       correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_upregulated <- GO_upregulated$result
# Save results.
write_xlsx(results_GO_upregulated, path = paste0("AMRvsIFTA_results_GO_upregulated.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
AMRvsIFTA_GO_up <- BP_GO_up_filtered
# Obtain top 3/10 GOs by lowest p_value.
GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
AMRvsIFTA_GO_top3  <- top_n(BP_GO_up_filtered , n = -5, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_up <- unique(GO_up_top10$term_name)
AMRvsIFTA_GO_top3_terms <- unique(AMRvsIFTA_GO_top3$term_name)



# Plot results
AMRvsIFTA_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                             y = factor(source), size = precision, fill = -log10(p_value),
                                             color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                             subtitle = "AMR vs IF/TA") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = AMRvsIFTA_up_plot, filename = "AMRvsIFTA_GO_up.jpeg")

### DOWNREGULATED ----
# Obtain downregulated gene names.
downregulated <- AMRvsIFTA$gene[AMRvsIFTA$avg_log2FC < 0 & 
                                  AMRvsIFTA$p_val_adj < 0.01]
# GO analysis using downregulated genes.
GO_downregulated <- gost(query = downregulated, 
                         organism = "hsapiens", ordered_query = T,
                         correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_downregulated <- GO_downregulated$result
# Save results.
write_xlsx(results_GO_downregulated, path = paste0("AMRvsIFTA_results_GO_downregulated.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_down_filtered <- filter(results_GO_downregulated, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_down_filtered <- dplyr::filter(GO_down_filtered, source == "GO:BP")
# Obtain top 10 GOs by lowest p_value.
GO_down_top10 <- top_n(BP_GO_down_filtered, n = -10, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_down <- unique(GO_down_top10$term_name)

# Plot results
AMRvsIFTA_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                 y = factor(source), size = precision, fill = -log10(p_value),
                                                 color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                             subtitle = "AMR vs IF/TA") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = AMRvsIFTA_down_plot, filename = "AMRvsIFTA_GO_down.jpeg")


### ALL ----
# Obtain all gene names.
all_genes <- AMRvsIFTA$gene[AMRvsIFTA$p_val_adj < 0.01]
# GO analysis using all genes.
GO_allgenes <- gost(query = all_genes, 
                    organism = "hsapiens", ordered_query = T,
                    correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_all <- GO_allgenes$result
# Save results.
write_xlsx(results_GO_all, path = paste0("AMRvsIFTA_results_GO_all.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_all_filtered <- filter(results_GO_all, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_all_filtered <- dplyr::filter(GO_all_filtered, source == "GO:BP")
# Obtain top 10 GOs by lowest p_value.
GO_all_top10 <- top_n(BP_GO_all_filtered, n = -10, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_all <- unique(GO_all_top10$term_name)

# Plot results
AMRvsIFTA_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                y = factor(source), size = precision, fill = -log10(p_value),
                                                color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                             subtitle = "AMR vs IF/TA") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = AMRvsIFTA_down_plot, filename = "AMRvsIFTA_GO_all.jpeg")

### Save results ----
AMRvsIFTA_GO_results <- NULL
AMRvsIFTA_GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
saveRDS(AMRvsIFTA_GO_results, file = paste0("gProfiler.RDS"))


# Plots top 10 of upregulated - downregulated - all
# UP
GO_upregulated$result <- BP_GO_up_filtered
plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                     subtitle = "AMR vs IF/TA")
publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_upregulated.jpeg"))

# DOWN
GO_downregulated$result <- BP_GO_down_filtered
plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                       subtitle = "AMR vs IF/TA")
publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_downregulated.jpeg"))
#ALL
GO_allgenes$result <- BP_GO_all_filtered
plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                  subtitle = "AMR vs IF/TA")
publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_all.jpeg"))
## Combined results ----
top3_terms <- c(TCMRvsIFTA_GO_top3_terms, AMRvsIFTA_GO_top3_terms)

TCMRvsIFTA_terms <- filter(TCMRvsIFTA_GO_up, term_name %in% top3_terms)
AMRvsIFTA_terms <- filter(AMRvsIFTA_GO_up, term_name %in% top3_terms)

TCMRvsIFTA_terms$diag <- "TCMR vs IF/TA"
AMRvsIFTA_terms$diag <- "AMR vs IF/TA"

test <- rbind(TCMRvsIFTA_terms, AMRvsIFTA_terms)

plot1 <- ggplot(test, aes(x = factor(term_name, level = top3_terms), 
                        y = factor(diag), size = precision, fill = -log10(p_value),
                        color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "   Top 5 GO:BP terms per disease comparison") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))

# Save
ggsave(plot = plot1, filename = "all_upregulated.jpeg", 
       path = "E:/MSc_EMC_project/Main_project/03_GO_analysis_outs/")

# STEP 5: GO Analysis per cluster. ----
## 5.1. TCMR vs IF/TA----
setwd(starting_directory)
# dir.create("percluster_TCMRvsIFTA")
setwd(paste0(starting_directory,"/percluster_TCMRvsIFTA"))

out_TCMRvsIFTA <- NULL

Idents(filtered_combined_samples) <- "clusterno_diag"
# Create output folder per cluster
for (cluster in levels(filtered_combined_samples@meta.data$seurat_clusters)) {
  dir.create(paste0(cluster))
  }

# GO analysis per cluster
for (cluster in levels(filtered_combined_samples@meta.data$seurat_clusters)) {
  # Set cluster no.
  cluster_no <- cluster
  print(paste0("Current cluster: ", cluster_no))
  if (cluster_no != 10){
    Idents(filtered_combined_samples) <- "clusterno_diag"
    TCMRvsIFTA <- FindMarkers(filtered_combined_samples,
                                     test.use = "MAST",
                                     min.pct = 0.25, logfc.threshold = 0.3,
                                     ident.1 = paste0(cluster_no, "_", TCMR), 
                                     ident.2 = paste0(cluster_no, "_", IFTA))
    
    # Sort by avg_log2FC.
    TCMRvsIFTA <- arrange(TCMRvsIFTA, -avg_log2FC)
    # Add genes as separate column.
    TCMRvsIFTA$gene <- rownames(TCMRvsIFTA)
    # Add cluster no.
    TCMRvsIFTA$cluster <- cluster_no
    
    # Create or add data to dataframe if cluster_no != 0.
    if (cluster_no == 0){
    out_TCMRvsIFTA <- TCMRvsIFTA
    } else {
      out_TCMRvsIFTA <- rbind(out_TCMRvsIFTA, TCMRvsIFTA)
    }
    
    # Save as Excel file.
    write_xlsx(TCMRvsIFTA, path = paste0(cluster_no, "/", "cl_TCMRvsIFTA.xlsx"))
    
    # Obtain top 10 top and top 10 bottom genes by Log2FC
    top_10_top_bottom <- c(rownames(out_TCMRvsIFTA)[1:10], # TOP
                           rownames(tail(out_TCMRvsIFTA, n = 10))) # BOTTOM
    Idents(filtered_combined_samples) <- "seurat_clusters"
    
    to_save <- VlnPlot(subset(filtered_combined_samples, idents = cluster_no), 
            features = top_10_top_bottom, group.by = "diag",
            pt.size = 0.1, col = palette.colors(palette = "Okabe-Ito")[4:8]) + RestoreLegend()
    ggsave(plot = to_save, paste0(cluster_no, "/", "cl_top_10_tb_VlnPlot.jpeg"),
           width = 25, height = 25)
    
    Idents(filtered_combined_samples) <- "clusterno_diag"
    
    ### UPREGULATED ----
    # Obtain upregulated gene names.
    upregulated <- TCMRvsIFTA$gene[TCMRvsIFTA$avg_log2FC > 0 & 
                                     TCMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using upregulated genes.
    GO_upregulated <- gost(query = upregulated, 
                           organism = "hsapiens", ordered_query = T,
                           correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_upregulated <- GO_upregulated$result
    # Save results.
    write_xlsx(results_GO_upregulated, path = paste0(cluster_no, "/","TCMRvsIFTA_results_GO_upregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
    assign(paste0("cl",cluster_no, "_TCMRvsIFTA_GO_up"), BP_GO_up_filtered)
    # Obtain top 10 GOs by lowest p_value.
    GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
    assign(paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5"), top_n(BP_GO_up_filtered, n = -5, wt = p_value))
    # Obtain top 10 term names.
    top10_term_names_GO_up <- unique(GO_up_top10$term_name)
    # paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5_terms") <- unique((paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5"))$term_name)
    assign(paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5_terms"), unique(get(paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5"))$term_name))
    
    # Plot results
    TCMRvsIFTA_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                                  y = factor(source), size = precision, fill = -log10(p_value),
                                                  color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                                 subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = TCMRvsIFTA_up_plot, filename = paste0(cluster_no, "/", "TCMRvsIFTA_GO_up.jpeg"))
    
    ### DOWNREGULATED ----
    # Obtain downregulated gene names.
    downregulated <- TCMRvsIFTA$gene[TCMRvsIFTA$avg_log2FC < 0 & 
                                       TCMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using downregulated genes.
    GO_downregulated <- gost(query = downregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_downregulated <- GO_downregulated$result
    # Save results.
    write_xlsx(results_GO_downregulated, path = paste0(cluster_no, "/", "TCMRvsIFTA_results_GO_downregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_down_filtered <- filter(results_GO_downregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_down_filtered <- dplyr::filter(GO_down_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_down_top10 <- top_n(BP_GO_down_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_down <- unique(GO_down_top10$term_name)
    
    # Plot results
    TCMRvsIFTA_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                      y = factor(source), size = precision, fill = -log10(p_value),
                                                      color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                                 subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = TCMRvsIFTA_down_plot, filename = paste0(cluster_no, "/", "TCMRvsIFTA_GO_down.jpeg"))
    
    ### ALL ----
    # Obtain all gene names.
    all_genes <- TCMRvsIFTA$gene[TCMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using all genes.
    GO_allgenes <- gost(query = all_genes, 
                        organism = "hsapiens", ordered_query = T,
                        correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_all <- GO_allgenes$result
    # Save results.
    write_xlsx(results_GO_all, path = paste0(cluster_no, "/", "TCMRvsIFTA_results_GO_all.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_all_filtered <- filter(results_GO_all, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_all_filtered <- dplyr::filter(GO_all_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_all_top10 <- top_n(BP_GO_all_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_all <- unique(GO_all_top10$term_name)
    
    # Plot results
    TCMRvsIFTA_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                     y = factor(source), size = precision, fill = -log10(p_value),
                                                     color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                                 subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = TCMRvsIFTA_down_plot, filename = paste0(cluster_no, "/", "TCMRvsIFTA_GO_all.jpeg"))
    
    ### Save results ----
    # Combine results
    if (cluster_no == 0){
      GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
    } else {
      GO_results_new <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
      GO_results <- rbind(GO_results, GO_results_new)
    }
    
    GO_upregulated$result <- BP_GO_up_filtered
    plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                         subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_upregulated.jpeg"))
    
    # DOWN
    GO_downregulated$result <- BP_GO_down_filtered
    plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                           subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_downregulated.jpeg"))
    #ALL
    GO_allgenes$result <- BP_GO_all_filtered
    plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                      subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_all.jpeg"))
    
    saveRDS(GO_results,
            file = paste0(cluster_no, "/", "cl_gProfiler.RDS"))
    saveRDS(out_TCMRvsIFTA,
            paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    write_xlsx(out_TCMRvsIFTA,
               paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    
    print(paste0("Cluster ", cluster_no, " done."))
    
  }
  else {
    print(paste0("!!! Cluster number is 10 !!!"))
    Idents(filtered_combined_samples) <- "clusterno_diag"
    TCMRvsIFTA <- FindMarkers(filtered_combined_samples,
                              test.use = "MAST",
                              min.pct = 0.25, logfc.threshold = 0.3,
                              ident.1 = paste0(cluster_no, "_", 'aTCMR2B'), 
                              ident.2 = paste0(cluster_no, "_", IFTA))
    
    # Sort by avg_log2FC.
    TCMRvsIFTA <- arrange(TCMRvsIFTA, -avg_log2FC)
    # Add genes as separate column.
    TCMRvsIFTA$gene <- rownames(TCMRvsIFTA)
    # Add cluster no.
    TCMRvsIFTA$cluster <- cluster_no
    
    # Create or add data to dataframe if cluster_no != 0.
    if (cluster_no == 0){
      out_TCMRvsIFTA <- TCMRvsIFTA
    } else {
      out_TCMRvsIFTA <- rbind(out_TCMRvsIFTA, TCMRvsIFTA)
    }
    
    # Save as Excel file.
    write_xlsx(TCMRvsIFTA, path = paste0(cluster_no, "/", "cl_TCMRvsIFTA.xlsx"))
    
    # Obtain top 10 top and top 10 bottom genes by Log2FC
    top_10_top_bottom <- c(rownames(out_TCMRvsIFTA)[1:10], # TOP
                           rownames(tail(out_TCMRvsIFTA, n = 10))) # BOTTOM
    Idents(filtered_combined_samples) <- "seurat_clusters"
    
    to_save <- VlnPlot(subset(filtered_combined_samples, idents = cluster_no), 
                       features = top_10_top_bottom, group.by = "diag",
                       pt.size = 0.1, col = palette.colors(palette = "Okabe-Ito")[4:8]) + RestoreLegend()
    ggsave(plot = to_save, paste0(cluster_no, "/", "cl_top_10_tb_VlnPlot.jpeg"),
           width = 25, height = 25)
    
    Idents(filtered_combined_samples) <- "clusterno_diag"
    
    ### UPREGULATED ----
    # Obtain upregulated gene names.
    upregulated <- TCMRvsIFTA$gene[TCMRvsIFTA$avg_log2FC > 0 & 
                                     TCMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using upregulated genes.
    GO_upregulated <- gost(query = upregulated, 
                           organism = "hsapiens", ordered_query = T,
                           correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_upregulated <- GO_upregulated$result
    # Save results.
    write_xlsx(results_GO_upregulated, path = paste0(cluster_no, "/","TCMRvsIFTA_results_GO_upregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
    assign(paste0("cl",cluster_no, "_TCMRvsIFTA_GO_up"), BP_GO_up_filtered)
    # Obtain top 10 GOs by lowest p_value.
    GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
    assign(paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5"), top_n(BP_GO_up_filtered, n = -5, wt = p_value))
    # Obtain top 10 term names.
    top10_term_names_GO_up <- unique(GO_up_top10$term_name)
    assign(paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5_terms"), unique(get(paste0("cl", cluster_no, "_TCMRvsIFTA_GO_top5"))$term_name))
    
    # Plot results
    TCMRvsIFTA_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                                  y = factor(source), size = precision, fill = -log10(p_value),
                                                  color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                                 subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = TCMRvsIFTA_up_plot, filename = paste0(cluster_no, "/", "TCMRvsIFTA_GO_up.jpeg"))
    
    
    # Combine results
    if (cluster_no == 0){
      GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
    } else {
      GO_results_new <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
      GO_results <- rbind(GO_results, GO_results_new)
    }
    
    ### DOWNREGULATED ----
    # Obtain downregulated gene names.
    downregulated <- TCMRvsIFTA$gene[TCMRvsIFTA$avg_log2FC < 0 & 
                                       TCMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using downregulated genes.
    GO_downregulated <- gost(query = downregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_downregulated <- GO_downregulated$result
    # Save results.
    write_xlsx(results_GO_downregulated, path = paste0(cluster_no, "/", "TCMRvsIFTA_results_GO_downregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_down_filtered <- filter(results_GO_downregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_down_filtered <- dplyr::filter(GO_down_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_down_top10 <- top_n(BP_GO_down_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_down <- unique(GO_down_top10$term_name)
    
    # Plot results
    TCMRvsIFTA_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                      y = factor(source), size = precision, fill = -log10(p_value),
                                                      color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                                 subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = TCMRvsIFTA_down_plot, filename = paste0(cluster_no, "/", "TCMRvsIFTA_GO_down.jpeg"))
    
    ### ALL ----
    # Obtain all gene names.
    all_genes <- TCMRvsIFTA$gene[TCMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using all genes.
    GO_allgenes <- gost(query = all_genes, 
                        organism = "hsapiens", ordered_query = T,
                        correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_all <- GO_allgenes$result
    # Save results.
    write_xlsx(results_GO_all, path = paste0(cluster_no, "/", "TCMRvsIFTA_results_GO_all.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_all_filtered <- filter(results_GO_all, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_all_filtered <- dplyr::filter(GO_all_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_all_top10 <- top_n(BP_GO_all_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_all <- unique(GO_all_top10$term_name)
    
    # Plot results
    TCMRvsIFTA_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                     y = factor(source), size = precision, fill = -log10(p_value),
                                                     color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                                 subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = TCMRvsIFTA_down_plot, filename = paste0(cluster_no, "/", "TCMRvsIFTA_GO_all.jpeg"))
    
    ### Save results ----
    # Combine results
    if (cluster_no == 0){
      GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
    } else {
      GO_results_new <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
      GO_results <- rbind(GO_results, GO_results_new)
    }
    
    GO_upregulated$result <- BP_GO_up_filtered
    plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                         subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_upregulated.jpeg"))
    
    # DOWN
    GO_downregulated$result <- BP_GO_down_filtered
    plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                           subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_downregulated.jpeg"))
    #ALL
    GO_allgenes$result <- BP_GO_all_filtered
    plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                      subtitle = paste0("TCMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_all.jpeg"))
    
    saveRDS(GO_results,
            file = paste0(cluster_no, "/", "cl_gProfiler.RDS"))
    saveRDS(out_TCMRvsIFTA,
            paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    write_xlsx(out_TCMRvsIFTA,
               paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    
    print(paste0("Cluster ", cluster_no, " done."))
  }
}
  
  
    
## Combined results ----
top5_terms_TCMRvsIFTA <- c()
for (cluster in levels(filtered_combined_samples@meta.data$seurat_clusters)){
  top5 <- get(paste0("cl", cluster, "_TCMRvsIFTA_GO_top5_terms"))
  top5_terms_TCMRvsIFTA <- c(top5_terms_TCMRvsIFTA, top5)
}
top5_terms_TCMRvsIFTA <- unique(top5_terms_TCMRvsIFTA)

df_TCMRvsIFTA <- data.frame()
for (cluster in levels(filtered_combined_samples@meta.data$seurat_clusters)){
  GO_up <- get(paste0("cl", cluster, "_TCMRvsIFTA_GO_up"))
  GO_up$cluster <- cluster
  df_TCMRvsIFTA <- rbind(df_TCMRvsIFTA, GO_up)
}

df_TCMRvsIFTA$diag <- "TCMR vs IF/TA"
TCMRvsIFTA_terms_percl <- filter(df_TCMRvsIFTA, term_name %in% top5_terms_TCMRvsIFTA)
TCMRvsIFTA_terms_percl$cluster <- as.numeric(TCMRvsIFTA_terms_percl$cluster)

plot <- ggplot(TCMRvsIFTA_terms_percl, aes(x = factor(term_name, level = top5_terms_TCMRvsIFTA), 
                         y = factor(cluster), size = precision, fill = -log10(p_value),
                         color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top 5 GO:BP terms per cluster", subtitle = "TCMR vs IF/TA") +
  theme(axis.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0, hjust = 0.4), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0))

ggsave(plot = plot, filename = "top5_per_cluster_TCMRvsIFTA.jpeg", 
       path = "E:/MSc_EMC_project/Main_project/03_GO_analysis_outs/",
       height = 10)


## 5.2. AMR vs IF/TA ----
setwd(starting_directory)
dir.create("percluster_AMRvsIFTA")
setwd(paste0(starting_directory,"/percluster_AMRvsIFTA"))

out_AMRvsIFTA <- NULL

for (cluster in levels(filtered_combined_samples@meta.data$seurat_clusters)) {
  dir.create(paste0(cluster))
}

#length of filtered_combined_samples$clusterno_diag[filtered_combined_samples$clusterno_diag == paste0(cluster_no, "_", AMR)]
# is 0, so cluster 10 is skipped

# GO analysis per cluster
for (cluster in levels(filtered_combined_samples@meta.data$seurat_clusters)) {
  # Set cluster no.
  cluster_no <- cluster
  print(paste0("Current cluster: ", cluster_no))
  if (cluster_no != 10){
    Idents(filtered_combined_samples) <- "clusterno_diag"
    AMRvsIFTA <- FindMarkers(filtered_combined_samples,
                             test.use = "MAST",
                             min.pct = 0.25, logfc.threshold = 0.3,
                             ident.1 = paste0(cluster_no, "_", AMR), 
                             ident.2 = paste0(cluster_no, "_", IFTA))
    
    # Sort by avg_log2FC.
    AMRvsIFTA <- arrange(AMRvsIFTA, -avg_log2FC)
    # Add genes as separate column.
    AMRvsIFTA$gene <- rownames(AMRvsIFTA)
    # Add cluster no.
    AMRvsIFTA$cluster <- cluster_no
    
    # Create or add data to dataframe if cluster_no != 0.
    if (cluster_no == 0){
      out_AMRvsIFTA <- AMRvsIFTA
    }
    else {
      out_AMRvsIFTA <- rbind(out_AMRvsIFTA, AMRvsIFTA)
    }
    
    # Save as Excel file.
    write_xlsx(AMRvsIFTA, path = paste0(cluster_no, "/", "cl_AMRvsIFTA.xlsx"))
    
    # Obtain top 10 top and top 10 bottom genes by Log2FC
    top_10_top_bottom <- c(rownames(out_AMRvsIFTA)[1:10], # TOP
                           rownames(tail(out_AMRvsIFTA, n = 10))) # BOTTOM
    Idents(filtered_combined_samples) <- "seurat_clusters"
    
    to_save <- VlnPlot(subset(filtered_combined_samples, idents = cluster_no), 
                       features = top_10_top_bottom, group.by = "diag",
                       pt.size = 0.1, col = palette.colors(palette = "Okabe-Ito")[4:8]) + RestoreLegend()
    ggsave(paste0(cluster_no, "/", "cl_top_10_tb_VlnPlot.jpeg"),
           width = 25, height = 25)
    
    Idents(filtered_combined_samples) <- "clusterno_diag"
    
    ### UPREGULATED ----
    # Obtain upregulated gene names.
    upregulated <- AMRvsIFTA$gene[AMRvsIFTA$avg_log2FC > 0 & 
                                    AMRvsIFTA$p_val_adj < 0.01]
    if (length(upregulated) > 5){
      # GO analysis using upregulated genes.
      GO_upregulated <- gost(query = upregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
      # Save results as dataframe.
      results_GO_upregulated <- GO_upregulated$result
      # Save results.
      write_xlsx(results_GO_upregulated, path = paste0(cluster_no, "/","AMRvsIFTA_results_GO_upregulated.xlsx"))
      
      # Plot
      # Filter out general GOs (term_size = no. of genes annotated to the term).
      GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
      # Only keep results from GO:BP source.
      BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
      assign(paste0("cl",cluster_no, "_AMRvsIFTA_GO_up"), BP_GO_up_filtered)
      # Obtain top 10 GOs by lowest p_value.
      GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
      assign(paste0("cl", cluster_no, "_AMRvsIFTA_GO_top5"), top_n(BP_GO_up_filtered, n = -5, wt = p_value))
      # Obtain top 10 term names.
      top10_term_names_GO_up <- unique(GO_up_top10$term_name)
      assign(paste0("cl", cluster_no, "_AMRvsIFTA_GO_top5_terms"), unique(get(paste0("cl", cluster_no, "_AMRvsIFTA_GO_top5"))$term_name))
      
      # Plot results
      AMRvsIFTA_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                                   y = factor(source), size = precision, fill = -log10(p_value),
                                                   color = -log10(p_value))) +
        geom_point() + 
        coord_flip() +
        ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                                   subtitle = paste0("AMR vs IF/TA, cluster ", cluster_no)) +
        theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
              legend.title = element_text(size = 8), plot.title.position = "plot", 
              plot.title = element_text(hjust = 0.0))
      # Save
      ggsave(plot = AMRvsIFTA_up_plot, filename = paste0(cluster_no, "/", "AMRvsIFTA_GO_up.jpeg"))
    }
    
    ### DOWNREGULATED ----
    # Obtain downregulated gene names.
    downregulated <- AMRvsIFTA$gene[AMRvsIFTA$avg_log2FC < 0 & 
                                      AMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using downregulated genes.
    GO_downregulated <- gost(query = downregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_downregulated <- GO_downregulated$result
    # Save results.
    write_xlsx(results_GO_downregulated, path = paste0(cluster_no, "/", "AMRvsIFTA_results_GO_downregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_down_filtered <- filter(results_GO_downregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_down_filtered <- dplyr::filter(GO_down_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_down_top10 <- top_n(BP_GO_down_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_down <- unique(GO_down_top10$term_name)
    
    # Plot results
    AMRvsIFTA_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                     y = factor(source), size = precision, fill = -log10(p_value),
                                                     color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                                 subtitle = paste0("AMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = AMRvsIFTA_down_plot, filename = paste0(cluster_no, "/", "AMRvsIFTA_GO_down.jpeg"))
    
    ### ALL ----
    # Obtain all gene names.
    all_genes <- AMRvsIFTA$gene[AMRvsIFTA$p_val_adj < 0.01]
    # GO analysis using all genes.
    GO_allgenes <- gost(query = all_genes, 
                        organism = "hsapiens", ordered_query = T,
                        correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_all <- GO_allgenes$result
    # Save results.
    write_xlsx(results_GO_all, path = paste0(cluster_no, "/", "AMRvsIFTA_results_GO_all.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_all_filtered <- filter(results_GO_all, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_all_filtered <- dplyr::filter(GO_all_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_all_top10 <- top_n(BP_GO_all_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_all <- unique(GO_all_top10$term_name)
    
    # Plot results
    AMRvsIFTA_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                    y = factor(source), size = precision, fill = -log10(p_value),
                                                    color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                                 subtitle = paste0("AMR vs IF/TA, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = AMRvsIFTA_down_plot, filename = paste0(cluster_no, "/", "AMRvsIFTA_GO_all.jpeg"))
    
    ### Save results ----
    # Combine results
    if (cluster_no == 0){
      GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
    } else {
      GO_results_new <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
      GO_results <- rbind(GO_results, GO_results_new)
    }
    
    GO_upregulated$result <- BP_GO_up_filtered
    plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                         subtitle = paste0("AMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_upregulated.jpeg"))
    
    # DOWN
    GO_downregulated$result <- BP_GO_down_filtered
    plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                           subtitle = paste0("AMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_downregulated.jpeg"))
    #ALL
    GO_allgenes$result <- BP_GO_all_filtered
    plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                      subtitle = paste0("AMR vs IF/TA, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_all.jpeg"))
    
    saveRDS(GO_results,
            file = paste0(cluster_no, "/", "cl_gProfiler.RDS"))
    saveRDS(out_AMRvsIFTA,
            paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    write_xlsx(out_AMRvsIFTA,
               paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    
    print(paste0("Cluster ", cluster_no, " done."))
  }
}

## Combined results ----
top5_terms_AMRvsIFTA <- c()
for (cluster in 1:9){
  top5 <- get(paste0("cl", cluster, "_AMRvsIFTA_GO_top5_terms"))
  top5_terms_AMRvsIFTA <- c(top5_terms_AMRvsIFTA, top5)
}
top5_terms_AMRvsIFTA <- unique(top5_terms_AMRvsIFTA)

df_AMRvsIFTA <- data.frame()
for (cluster in 1:9){
  GO_up <- get(paste0("cl", cluster, "_AMRvsIFTA_GO_up"))
  GO_up$cluster <- cluster
  df_AMRvsIFTA <- rbind(df_AMRvsIFTA, GO_up)
}

df_AMRvsIFTA$diag <- "AMR vs IF/TA"
AMRvsIFTA_terms_percl <- filter(df_AMRvsIFTA, term_name %in% top5_terms_AMRvsIFTA)
AMRvsIFTA_terms_percl$cluster <- as.numeric(AMRvsIFTA_terms_percl$cluster)


plot <- ggplot(AMRvsIFTA_terms_percl, aes(x = factor(term_name, level = top5_terms_AMRvsIFTA), 
                                           y = factor(cluster), size = precision, fill = -log10(p_value),
                                           color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top 5 GO:BP terms per cluster", subtitle = "AMR vs IF/TA") +
  theme(axis.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0, hjust = 0.4), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0))

ggsave(plot = plot, filename = "top5_per_cluster_AMRvsIFTA.jpeg", 
       path = "E:/MSc_EMC_project/Main_project/03_GO_analysis_outs/",
       height = 10)

## Both diseases combined ----
allterms <- c(top5_terms_TCMRvsIFTA, top5_terms_AMRvsIFTA)
all <- rbind(TCMRvsIFTA_terms_percl, AMRvsIFTA_terms_percl)


plot2 <- ggplot(all, aes(x = factor(term_name, level = allterms), 
                                          y = factor(cluster), size = precision, fill = -log10(p_value),
                                          color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  facet_wrap(~diag) +
  ylab("") + xlab("") + labs(title = "   Top 5 GO:BP terms per cluster in AMR and TCMR vs IF/TA") +
  theme(axis.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0, hjust = 0.4), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0))

ggsave(plot = plot2, filename = "top5_per_cluster_alldiseases.jpeg", 
       path = "E:/MSc_EMC_project/Main_project/03_GO_analysis_outs/",
       width = 10)

## Final figure ----
Fig_5 <- ggarrange(plot1, plot2, labels = c("A", "B"), widths = c(1,2))
ggsave(plot = Fig_5, filename = "GO_results.jpeg", path = "E:/MSc_EMC_project/Main_project/03_GO_analysis_outs/", width = 15)
