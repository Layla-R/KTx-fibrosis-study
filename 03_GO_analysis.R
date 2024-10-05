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
## 4.1. IFTA vs TCMR ----
setwd(starting_directory)
# dir.create("IFTAvsTCMR")
setwd(paste0(starting_directory, "/IFTAvsTCMR"))

Idents(filtered_combined_samples) <- 'diag'
IFTAvsTCMR <- FindMarkers(filtered_combined_samples,
                          test.use = "MAST",
                          min.pct = 0.25, logfc.threshold = 0.3,
                          ident.1 = IFTA, 
                          ident.2 = TCMR)

# Sort by avg_log2FC.
IFTAvsTCMR <- arrange(IFTAvsTCMR, -avg_log2FC)
# Add gene names as separate column.
IFTAvsTCMR$gene <- rownames(IFTAvsTCMR)
# Save as Excel file & RDS.
write_xlsx(IFTAvsTCMR, path = paste0("DE_FindMarkers_gProfiler_IFTAvsTCMR.xlsx"))
saveRDS(IFTAvsTCMR, paste0("DE_FindMarkers_gProfiler.RDS"))

### UPREGULATED ----
# Obtain upregulated gene names.
upregulated <- IFTAvsTCMR$gene[IFTAvsTCMR$avg_log2FC > 0 & 
                                 IFTAvsTCMR$p_val_adj < 0.01]
# GO analysis using upregulated genes.
GO_upregulated <- gost(query = upregulated, 
                       organism = "hsapiens", ordered_query = T,
                       correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_upregulated <- GO_upregulated$result
# Save results.
write_xlsx(results_GO_upregulated, path = paste0("IFTAvsTCMR_results_GO_upregulated.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
# Obtain top 10 GOs by lowest p_value.
GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_up <- unique(GO_up_top10$term_name)

# Plot results
IFTAvsTCMR_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                              y = factor(source), size = precision, fill = -log10(p_value),
                                              color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                             subtitle = "IF/TA vs TCMR") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = IFTAvsTCMR_up_plot, filename = "IFTAvsTCMR_GO_up.jpeg")


### DOWNREGULATED ----
# Obtain downregulated gene names.
downregulated <- IFTAvsTCMR$gene[IFTAvsTCMR$avg_log2FC < 0 & 
                                   IFTAvsTCMR$p_val_adj < 0.01]
# GO analysis using downregulated genes.
GO_downregulated <- gost(query = downregulated, 
                         organism = "hsapiens", ordered_query = T,
                         correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_downregulated <- GO_downregulated$result
# Save results.
write_xlsx(results_GO_downregulated, path = paste0("IFTAvsTCMR_results_GO_downregulated.xlsx"))

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
IFTAvsTCMR_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                  y = factor(source), size = precision, fill = -log10(p_value),
                                                  color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                             subtitle = "IF/TA vs TCMR") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = IFTAvsTCMR_down_plot, filename = "IFTAvsTCMR_GO_down.jpeg")

### ALL ----
# Obtain all gene names.
all_genes <- IFTAvsTCMR$gene[IFTAvsTCMR$p_val_adj < 0.01]
# GO analysis using all genes.
GO_allgenes <- gost(query = all_genes, 
                    organism = "hsapiens", ordered_query = T,
                    correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_all <- GO_allgenes$result
# Save results.
write_xlsx(results_GO_all, path = paste0("IFTAvsTCMR_results_GO_all.xlsx"))

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
IFTAvsTCMR_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                 y = factor(source), size = precision, fill = -log10(p_value),
                                                 color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                             subtitle = "IF/TA vs TCMR") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = IFTAvsTCMR_down_plot, filename = "IFTAvsTCMR_GO_all.jpeg")

### Save results ----
IFTAvsTCMR_GO_results <- NULL
IFTAvsTCMR_GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
saveRDS(IFTAvsTCMR_GO_results,
        file = paste0("gProfiler.RDS"))

# Plots top 10 of upregulated - downregulated - all
# UP
GO_upregulated$result <- BP_GO_up_filtered
plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                     subtitle = "IF/TA vs TCMR")
publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_upregulated.jpeg"))

# DOWN
GO_downregulated$result <- BP_GO_down_filtered
plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                       subtitle = "IF/TA vs TCMR")
publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_downregulated.jpeg"))
#ALL
GO_allgenes$result <- BP_GO_all_filtered
plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                  subtitle = "IF/TA vs TCMR")
publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_all.jpeg"))

## 4.2. IFTA vs AMR ----
setwd(starting_directory)
# dir.create("IFTAvsAMR")
setwd(paste0(starting_directory, "/IFTAvsAMR"))

Idents(filtered_combined_samples) <- 'diag'
IFTAvsAMR <- FindMarkers(filtered_combined_samples,
                         test.use = "MAST",
                         min.pct = 0.25, logfc.threshold = 0.3,
                         ident.1 = IFTA, 
                         ident.2 = AMR)

# Sort by avg_log2FC.
IFTAvsAMR <- arrange(IFTAvsTCMR, -avg_log2FC)
# Add genes as separate column.
IFTAvsAMR$gene <- rownames(IFTAvsTCMR)
# Save as Excel file & RDS.
write_xlsx(IFTAvsAMR, path = paste0("DE_FindMarkers_gProfiler_IFTAvsAMR.xlsx"))
saveRDS(IFTAvsAMR, paste0("DE_FindMarkers_gProfiler.RDS"))

### UPREGULATED ----
# Obtain upregulated gene names.
upregulated <- IFTAvsAMR$gene[IFTAvsAMR$avg_log2FC > 0 & 
                                IFTAvsAMR$p_val_adj < 0.01]
# GO analysis using upregulated genes.
GO_upregulated <- gost(query = upregulated, 
                       organism = "hsapiens", ordered_query = T,
                       correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_upregulated <- GO_upregulated$result
# Save results.
write_xlsx(results_GO_upregulated, path = paste0("IFTAvsAMR_results_GO_upregulated.xlsx"))

# Plot
# Filter out general GOs (term_size = no. of genes annotated to the term).
GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
# Only keep results from GO:BP source.
BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
# Obtain top 10 GOs by lowest p_value.
GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
# Obtain top 10 term names.
top10_term_names_GO_up <- unique(GO_up_top10$term_name)

# Plot results
IFTAvsAMR_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                             y = factor(source), size = precision, fill = -log10(p_value),
                                             color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                             subtitle = "IF/TA vs AMR") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = IFTAvsAMR_up_plot, filename = "IFTAvsAMR_GO_up.jpeg")

### DOWNREGULATED ----
# Obtain downregulated gene names.
downregulated <- IFTAvsAMR$gene[IFTAvsAMR$avg_log2FC < 0 & 
                                  IFTAvsAMR$p_val_adj < 0.01]
# GO analysis using downregulated genes.
GO_downregulated <- gost(query = downregulated, 
                         organism = "hsapiens", ordered_query = T,
                         correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_downregulated <- GO_downregulated$result
# Save results.
write_xlsx(results_GO_downregulated, path = paste0("IFTAvsAMR_results_GO_downregulated.xlsx"))

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
IFTAvsAMR_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                 y = factor(source), size = precision, fill = -log10(p_value),
                                                 color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                             subtitle = "IF/TA vs AMR") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = IFTAvsAMR_down_plot, filename = "IFTAvsAMR_GO_down.jpeg")


### ALL ----
# Obtain all gene names.
all_genes <- IFTAvsAMR$gene[IFTAvsAMR$p_val_adj < 0.01]
# GO analysis using all genes.
GO_allgenes <- gost(query = all_genes, 
                    organism = "hsapiens", ordered_query = T,
                    correction_method = "g_SCS", domain_scope = "annotated")
# Save results as dataframe.
results_GO_all <- GO_allgenes$result
# Save results.
write_xlsx(results_GO_all, path = paste0("IFTAvsAMR_results_GO_all.xlsx"))

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
IFTAvsAMR_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                y = factor(source), size = precision, fill = -log10(p_value),
                                                color = -log10(p_value))) +
  geom_point() + 
  coord_flip() +
  ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                             subtitle = "IF/TA vs AMR") +
  theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 8), plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.0))
# Save
ggsave(plot = IFTAvsAMR_down_plot, filename = "IFTAvsAMR_GO_all.jpeg")

### Save results ----
IFTAvsAMR_GO_results <- NULL
IFTAvsAMR_GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
saveRDS(IFTAvsAMR_GO_results, file = paste0("gProfiler.RDS"))


# Plots top 10 of upregulated - downregulated - all
# UP
GO_upregulated$result <- BP_GO_up_filtered
plot <- gostplot(GO_upregulated, capped = F, interactive = F) + labs(title = "GSE of upregulated genes (Top 10)", 
                                                                     subtitle = "IF/TA vs AMR")
publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_upregulated.jpeg"))

# DOWN
GO_downregulated$result <- BP_GO_down_filtered
plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                       subtitle = "IF/TA vs AMR")
publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_downregulated.jpeg"))
#ALL
GO_allgenes$result <- BP_GO_all_filtered
plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                  subtitle = "IF/TA vs AMR")
publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                 width = 15, height = 12,
                 filename = paste0("GO_top10_all.jpeg"))
# STEP 5: GO Analysis per cluster. ----
## 5.1. IFTA vs TCMR ----
setwd(starting_directory)
# dir.create("percluster_IFTAvsTCMR")
setwd(paste0(starting_directory,"/percluster_IFTAvsTCMR"))

out_IFTAvsTCMR <- NULL

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
    IFTAvsTCMR <- FindMarkers(filtered_combined_samples,
                                     test.use = "MAST",
                                     min.pct = 0.25, logfc.threshold = 0.3,
                                     ident.1 = paste0(cluster_no, "_", IFTA), 
                                     ident.2 = paste0(cluster_no, "_", TCMR))
    
    # Sort by avg_log2FC.
    IFTAvsTCMR <- arrange(IFTAvsTCMR, -avg_log2FC)
    # Add genes as separate column.
    IFTAvsTCMR$gene <- rownames(IFTAvsTCMR)
    # Add cluster no.
    IFTAvsTCMR$cluster <- cluster_no
    
    # Create or add data to dataframe if cluster_no != 0.
    if (cluster_no == 0){
    out_IFTAvsTCMR <- IFTAvsTCMR
    } else {
      out_IFTAvsTCMR <- rbind(out_IFTAvsTCMR, IFTAvsTCMR)
    }
    
    # Save as Excel file.
    write_xlsx(IFTAvsTCMR, path = paste0(cluster_no, "/", "cl_IFTAvsTCMR.xlsx"))
    
    # Obtain top 10 top and top 10 bottom genes by Log2FC
    top_10_top_bottom <- c(rownames(out_IFTAvsTCMR)[1:10], # TOP
                           rownames(tail(out_IFTAvsTCMR, n = 10))) # BOTTOM
    Idents(filtered_combined_samples) <- "seurat_clusters"
    
    to_save <- VlnPlot(subset(filtered_combined_samples, idents = cluster_no), 
            features = top_10_top_bottom, group.by = "diag",
            pt.size = 0.1, col = palette.colors(palette = "Okabe-Ito")[4:8]) + RestoreLegend()
    ggsave(plot = to_save, paste0(cluster_no, "/", "cl_top_10_tb_VlnPlot.jpeg"),
           width = 25, height = 25)
    
    Idents(filtered_combined_samples) <- "clusterno_diag"
    
    ### UPREGULATED ----
    # Obtain upregulated gene names.
    upregulated <- IFTAvsTCMR$gene[IFTAvsTCMR$avg_log2FC > 0 & 
                                     IFTAvsTCMR$p_val_adj < 0.01]
    # GO analysis using upregulated genes.
    GO_upregulated <- gost(query = upregulated, 
                           organism = "hsapiens", ordered_query = T,
                           correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_upregulated <- GO_upregulated$result
    # Save results.
    write_xlsx(results_GO_upregulated, path = paste0(cluster_no, "/","IFTAvsTCMR_results_GO_upregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_up <- unique(GO_up_top10$term_name)
    
    # Plot results
    IFTAvsTCMR_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                                  y = factor(source), size = precision, fill = -log10(p_value),
                                                  color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                                 subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsTCMR_up_plot, filename = paste0(cluster_no, "/", "IFTAvsTCMR_GO_up.jpeg"))
    
    ### DOWNREGULATED ----
    # Obtain downregulated gene names.
    downregulated <- IFTAvsTCMR$gene[IFTAvsTCMR$avg_log2FC < 0 & 
                                       IFTAvsTCMR$p_val_adj < 0.01]
    # GO analysis using downregulated genes.
    GO_downregulated <- gost(query = downregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_downregulated <- GO_downregulated$result
    # Save results.
    write_xlsx(results_GO_downregulated, path = paste0(cluster_no, "/", "IFTAvsTCMR_results_GO_downregulated.xlsx"))
    
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
    IFTAvsTCMR_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                      y = factor(source), size = precision, fill = -log10(p_value),
                                                      color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                                 subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsTCMR_down_plot, filename = paste0(cluster_no, "/", "IFTAvsTCMR_GO_down.jpeg"))
    
    ### ALL ----
    # Obtain all gene names.
    all_genes <- IFTAvsTCMR$gene[IFTAvsTCMR$p_val_adj < 0.01]
    # GO analysis using all genes.
    GO_allgenes <- gost(query = all_genes, 
                        organism = "hsapiens", ordered_query = T,
                        correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_all <- GO_allgenes$result
    # Save results.
    write_xlsx(results_GO_all, path = paste0(cluster_no, "/", "IFTAvsTCMR_results_GO_all.xlsx"))
    
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
    IFTAvsTCMR_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                     y = factor(source), size = precision, fill = -log10(p_value),
                                                     color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                                 subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsTCMR_down_plot, filename = paste0(cluster_no, "/", "IFTAvsTCMR_GO_all.jpeg"))
    
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
                                                                         subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_upregulated.jpeg"))
    
    # DOWN
    GO_downregulated$result <- BP_GO_down_filtered
    plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                           subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_downregulated.jpeg"))
    #ALL
    GO_allgenes$result <- BP_GO_all_filtered
    plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                      subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_all.jpeg"))
    
    saveRDS(GO_results,
            file = paste0(cluster_no, "/", "cl_gProfiler.RDS"))
    saveRDS(out_IFTAvsTCMR,
            paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    write_xlsx(out_IFTAvsTCMR,
               paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    
    print(paste0("Cluster ", cluster_no, " done."))
    
  }
  else {
    print(paste0("!!! Cluster number is 10 !!!"))
    Idents(filtered_combined_samples) <- "clusterno_diag"
    IFTAvsTCMR <- FindMarkers(filtered_combined_samples,
                              test.use = "MAST",
                              min.pct = 0.25, logfc.threshold = 0.3,
                              ident.1 = paste0(cluster_no, "_", IFTA), 
                              ident.2 = paste0(cluster_no, "_", 'aTCMR2B'))
    
    # Sort by avg_log2FC.
    IFTAvsTCMR <- arrange(IFTAvsTCMR, -avg_log2FC)
    # Add genes as separate column.
    IFTAvsTCMR$gene <- rownames(IFTAvsTCMR)
    # Add cluster no.
    IFTAvsTCMR$cluster <- cluster_no
    
    # Create or add data to dataframe if cluster_no != 0.
    if (cluster_no == 0){
      out_IFTAvsTCMR <- IFTAvsTCMR
    } else {
      out_IFTAvsTCMR <- rbind(out_IFTAvsTCMR, IFTAvsTCMR)
    }
    
    # Save as Excel file.
    write_xlsx(IFTAvsTCMR, path = paste0(cluster_no, "/", "cl_IFTAvsTCMR.xlsx"))
    
    # Obtain top 10 top and top 10 bottom genes by Log2FC
    top_10_top_bottom <- c(rownames(out_IFTAvsTCMR)[1:10], # TOP
                           rownames(tail(out_IFTAvsTCMR, n = 10))) # BOTTOM
    Idents(filtered_combined_samples) <- "seurat_clusters"
    
    to_save <- VlnPlot(subset(filtered_combined_samples, idents = cluster_no), 
                       features = top_10_top_bottom, group.by = "diag",
                       pt.size = 0.1, col = palette.colors(palette = "Okabe-Ito")[4:8]) + RestoreLegend()
    ggsave(plot = to_save, paste0(cluster_no, "/", "cl_top_10_tb_VlnPlot.jpeg"),
           width = 25, height = 25)
    
    Idents(filtered_combined_samples) <- "clusterno_diag"
    
    ### UPREGULATED ----
    # Obtain upregulated gene names.
    upregulated <- IFTAvsTCMR$gene[IFTAvsTCMR$avg_log2FC > 0 & 
                                     IFTAvsTCMR$p_val_adj < 0.01]
    # GO analysis using upregulated genes.
    GO_upregulated <- gost(query = upregulated, 
                           organism = "hsapiens", ordered_query = T,
                           correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_upregulated <- GO_upregulated$result
    # Save results.
    write_xlsx(results_GO_upregulated, path = paste0(cluster_no, "/","IFTAvsTCMR_results_GO_upregulated.xlsx"))
    
    # Plot
    # Filter out general GOs (term_size = no. of genes annotated to the term).
    GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
    # Only keep results from GO:BP source.
    BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
    # Obtain top 10 GOs by lowest p_value.
    GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
    # Obtain top 10 term names.
    top10_term_names_GO_up <- unique(GO_up_top10$term_name)
    
    # Plot results
    IFTAvsTCMR_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                                  y = factor(source), size = precision, fill = -log10(p_value),
                                                  color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                                 subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsTCMR_up_plot, filename = paste0(cluster_no, "/", "IFTAvsTCMR_GO_up.jpeg"))
    
    
    # Combine results
    if (cluster_no == 0){
      GO_results <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
    } else {
      GO_results_new <- rbind(results_GO_upregulated, results_GO_downregulated, results_GO_all)
      GO_results <- rbind(GO_results, GO_results_new)
    }
    
    ### DOWNREGULATED ----
    # Obtain downregulated gene names.
    downregulated <- IFTAvsTCMR$gene[IFTAvsTCMR$avg_log2FC < 0 & 
                                       IFTAvsTCMR$p_val_adj < 0.01]
    # GO analysis using downregulated genes.
    GO_downregulated <- gost(query = downregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_downregulated <- GO_downregulated$result
    # Save results.
    write_xlsx(results_GO_downregulated, path = paste0(cluster_no, "/", "IFTAvsTCMR_results_GO_downregulated.xlsx"))
    
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
    IFTAvsTCMR_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                      y = factor(source), size = precision, fill = -log10(p_value),
                                                      color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                                 subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsTCMR_down_plot, filename = paste0(cluster_no, "/", "IFTAvsTCMR_GO_down.jpeg"))
    
    ### ALL ----
    # Obtain all gene names.
    all_genes <- IFTAvsTCMR$gene[IFTAvsTCMR$p_val_adj < 0.01]
    # GO analysis using all genes.
    GO_allgenes <- gost(query = all_genes, 
                        organism = "hsapiens", ordered_query = T,
                        correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_all <- GO_allgenes$result
    # Save results.
    write_xlsx(results_GO_all, path = paste0(cluster_no, "/", "IFTAvsTCMR_results_GO_all.xlsx"))
    
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
    IFTAvsTCMR_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                     y = factor(source), size = precision, fill = -log10(p_value),
                                                     color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                                 subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsTCMR_down_plot, filename = paste0(cluster_no, "/", "IFTAvsTCMR_GO_all.jpeg"))
    
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
                                                                         subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_upregulated.jpeg"))
    
    # DOWN
    GO_downregulated$result <- BP_GO_down_filtered
    plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                           subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_downregulated.jpeg"))
    #ALL
    GO_allgenes$result <- BP_GO_all_filtered
    plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                      subtitle = paste0("IF/TA vs TCMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_all.jpeg"))
    
    saveRDS(GO_results,
            file = paste0(cluster_no, "/", "cl_gProfiler.RDS"))
    saveRDS(out_IFTAvsTCMR,
            paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    write_xlsx(out_IFTAvsTCMR,
               paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    
    print(paste0("Cluster ", cluster_no, " done."))
  }
}
  
  
    
## 5.2. IFTA vs AMR ----
setwd(starting_directory)
# dir.create("percluster_IFTAvsAMR")
setwd(paste0(starting_directory,"/percluster_IFTAvsAMR"))

out_IFTAvsAMR <- NULL

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
    IFTAvsAMR <- FindMarkers(filtered_combined_samples,
                             test.use = "MAST",
                             min.pct = 0.25, logfc.threshold = 0.3,
                             ident.1 = paste0(cluster_no, "_", IFTA), 
                             ident.2 = paste0(cluster_no, "_", AMR))
    
    # Sort by avg_log2FC.
    IFTAvsAMR <- arrange(IFTAvsAMR, -avg_log2FC)
    # Add genes as separate column.
    IFTAvsAMR$gene <- rownames(IFTAvsAMR)
    # Add cluster no.
    IFTAvsAMR$cluster <- cluster_no
    
    # Create or add data to dataframe if cluster_no != 0.
    if (cluster_no == 0){
      out_IFTAvsAMR <- IFTAvsAMR
    }
    else {
      out_IFTAvsAMR <- rbind(out_IFTAvsAMR, IFTAvsAMR)
    }
    
    # Save as Excel file.
    write_xlsx(IFTAvsAMR, path = paste0(cluster_no, "/", "cl_IFTAvsAMR.xlsx"))
    
    # Obtain top 10 top and top 10 bottom genes by Log2FC
    top_10_top_bottom <- c(rownames(out_IFTAvsAMR)[1:10], # TOP
                           rownames(tail(out_IFTAvsAMR, n = 10))) # BOTTOM
    Idents(filtered_combined_samples) <- "seurat_clusters"
    
    to_save <- VlnPlot(subset(filtered_combined_samples, idents = cluster_no), 
                       features = top_10_top_bottom, group.by = "diag",
                       pt.size = 0.1, col = palette.colors(palette = "Okabe-Ito")[4:8]) + RestoreLegend()
    ggsave(paste0(cluster_no, "/", "cl_top_10_tb_VlnPlot.jpeg"),
           width = 25, height = 25)
    
    Idents(filtered_combined_samples) <- "clusterno_diag"
    
    ### UPREGULATED ----
    # Obtain upregulated gene names.
    upregulated <- IFTAvsAMR$gene[IFTAvsAMR$avg_log2FC > 0 & 
                                    IFTAvsAMR$p_val_adj < 0.01]
    if (length(upregulated) > 5){
      # GO analysis using upregulated genes.
      GO_upregulated <- gost(query = upregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
      # Save results as dataframe.
      results_GO_upregulated <- GO_upregulated$result
      # Save results.
      write_xlsx(results_GO_upregulated, path = paste0(cluster_no, "/","IFTAvsAMR_results_GO_upregulated.xlsx"))
      
      # Plot
      # Filter out general GOs (term_size = no. of genes annotated to the term).
      GO_up_filtered <- filter(results_GO_upregulated, term_size <= 2000)
      # Only keep results from GO:BP source.
      BP_GO_up_filtered <- dplyr::filter(GO_up_filtered, source == "GO:BP")
      # Obtain top 10 GOs by lowest p_value.
      GO_up_top10 <- top_n(BP_GO_up_filtered, n = -10, wt = p_value)
      # Obtain top 10 term names.
      top10_term_names_GO_up <- unique(GO_up_top10$term_name)
      
      # Plot results
      IFTAvsAMR_up_plot <- ggplot(GO_up_top10, aes(x = factor(term_name, level = top10_term_names_GO_up), 
                                                   y = factor(source), size = precision, fill = -log10(p_value),
                                                   color = -log10(p_value))) +
        geom_point() + 
        coord_flip() +
        ylab("") + xlab("") + labs(title = "Top GO:BP based on upregulated DE Proteins",
                                   subtitle = paste0("IF/TA vs AMR, cluster ", cluster_no)) +
        theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
              legend.title = element_text(size = 8), plot.title.position = "plot", 
              plot.title = element_text(hjust = 0.0))
      # Save
      ggsave(plot = IFTAvsAMR_up_plot, filename = paste0(cluster_no, "/", "IFTAvsAMR_GO_up.jpeg"))
    }
    
    ### DOWNREGULATED ----
    # Obtain downregulated gene names.
    downregulated <- IFTAvsAMR$gene[IFTAvsAMR$avg_log2FC < 0 & 
                                      IFTAvsAMR$p_val_adj < 0.01]
    # GO analysis using downregulated genes.
    GO_downregulated <- gost(query = downregulated, 
                             organism = "hsapiens", ordered_query = T,
                             correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_downregulated <- GO_downregulated$result
    # Save results.
    write_xlsx(results_GO_downregulated, path = paste0(cluster_no, "/", "IFTAvsAMR_results_GO_downregulated.xlsx"))
    
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
    IFTAvsAMR_down_plot <- ggplot(GO_down_top10, aes(x = factor(term_name, level = top10_term_names_GO_down), 
                                                     y = factor(source), size = precision, fill = -log10(p_value),
                                                     color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on downregulated DE Proteins",
                                 subtitle = paste0("IF/TA vs AMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsAMR_down_plot, filename = paste0(cluster_no, "/", "IFTAvsAMR_GO_down.jpeg"))
    
    ### ALL ----
    # Obtain all gene names.
    all_genes <- IFTAvsAMR$gene[IFTAvsAMR$p_val_adj < 0.01]
    # GO analysis using all genes.
    GO_allgenes <- gost(query = all_genes, 
                        organism = "hsapiens", ordered_query = T,
                        correction_method = "g_SCS", domain_scope = "annotated")
    # Save results as dataframe.
    results_GO_all <- GO_allgenes$result
    # Save results.
    write_xlsx(results_GO_all, path = paste0(cluster_no, "/", "IFTAvsAMR_results_GO_all.xlsx"))
    
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
    IFTAvsAMR_down_plot <- ggplot(GO_all_top10, aes(x = factor(term_name, level = top10_term_names_GO_all), 
                                                    y = factor(source), size = precision, fill = -log10(p_value),
                                                    color = -log10(p_value))) +
      geom_point() + 
      coord_flip() +
      ylab("") + xlab("") + labs(title = "Top GO:BP based on all DE Proteins",
                                 subtitle = paste0("IF/TA vs AMR, cluster ", cluster_no)) +
      theme(axis.text.y = , axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.title = element_text(size = 8), plot.title.position = "plot", 
            plot.title = element_text(hjust = 0.0))
    # Save
    ggsave(plot = IFTAvsAMR_down_plot, filename = paste0(cluster_no, "/", "IFTAvsTAMR_GO_all.jpeg"))
    
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
                                                                         subtitle = paste0("IF/TA vs AMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_upregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_upregulated.jpeg"))
    
    # DOWN
    GO_downregulated$result <- BP_GO_down_filtered
    plot <- gostplot(GO_downregulated, capped = F, interactive = F) + labs(title = "GSE of downregulated genes (Top 10)", 
                                                                           subtitle = paste0("IF/TA vs AMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_downregulated$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_downregulated.jpeg"))
    #ALL
    GO_allgenes$result <- BP_GO_all_filtered
    plot <- gostplot(GO_allgenes, capped = F, interactive = F) + labs(title = "GSE of all genes (Top 10)", 
                                                                      subtitle = paste0("IF/TA vs AMR, cluster ", cluster_no))
    publish_gostplot(plot, highlight_terms = GO_allgenes$result$term_id[1:10],
                     width = 15, height = 12,
                     filename = paste0(cluster_no, "/","GO_top10_all.jpeg"))
    
    saveRDS(GO_results,
            file = paste0(cluster_no, "/", "cl_gProfiler.RDS"))
    saveRDS(out_IFTAvsAMR,
            paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    write_xlsx(out_IFTAvsAMR,
               paste0(cluster_no, "/", "cl_DE_FindMarkers_gProfiler.RDS"))
    
    print(paste0("Cluster ", cluster_no, " done."))
  }
}
