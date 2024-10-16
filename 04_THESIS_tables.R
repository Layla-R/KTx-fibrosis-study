# Libraries ----
library(dplyr) 
library(gt)
library(gtsummary)
library(readxl)
library(tidyr)

# Table 1 ----
# Set WD
setwd("E:/MSc_EMC_project/Main_project/01_QC_outs/")

info_samples.xlsx <- read_excel("E:/MSc_EMC_project/Main_project/Main File_snRNAseq cases_v3xlsx.xlsx",
                                sheet = "Sheet2")
info_samples.xlsx <- info_samples.xlsx[, c(1, 2, 21:31)] 

colnames(info_samples.xlsx) <- c("PA_number", "Diagnosis", "ZIS", "Gender")
colnames(info_samples.xlsx)[13] <- "Gender"

sample_list <- c("H16-8962", "H14-11604", "H14-23350", 
                 "H16-9676", "H14-20579", "H14-6512", 
                 "H16-9416", "H14-6925", "H16-8593", 
                 "H13-22010", "H16-13550", "H16-5999", "H14-11792")
# Get rows only with used samples
sample_info <- info_samples.xlsx %>% filter(info_samples.xlsx[[1]] %in% sample_list)

# Adjustments
sample_info$Age <- as.numeric(as.character(sample_info$Age))
sample_info[sample_info == "aAMR, C4d+"] <- "AMR"
sample_info[sample_info == "aTCMR1B"] <- "TCMR"
sample_info[sample_info == "aTCMR2B"] <- "TCMR"
sample_info[sample_info == "IFTA"] <- "IF/TA"

# Create table
study_population <- sample_info |>
  select(Diagnosis, Age, Gender) |>
  tbl_summary(by = Diagnosis,
              type = all_continuous() ~ "continuous",
              statistic = list(
                all_categorical() ~ "{n}/{N} ({p}%)",
              all_continuous() ~ "{mean} ({sd})")) |>
  add_p()

# To save
study_population <- study_population %>%
  as_gt() %>%
  # Adjust options
  gt::tab_options(table.font.names = "Cambria", table.font.size = 10) %>%
  # Save
  gt::gtsave(filename = "Table1_study_population.png")


# Table 2 ----
setwd("E:/MSc_EMC_project/Main_project")
Visium_filenames <- c("BIO1_KH119_A", "BIO2_KH120_B", "BIO3_KH121_A",
                      "BIO4_KH122_B", "BIO5_KH123_A", "BIO6_KH124_B")
# Collect outputs
SR_output <- list()
for (filename in Visium_filenames){
  SR_output[[`filename`]] <- read.csv(paste0("//wsl.localhost/Ubuntu/home/layla7x/MSc_project/outs/",
                                             filename, "/outs/metrics_summary.csv"))
}

# Create data frame
SR_df <- NULL
for (i in 1:length(SR_output)){
  SR_df <- rbind(SR_df, SR_output[[i]])
}

# Adjustments
colnames(SR_df) <- gsub("\\.", " ", colnames(SR_df))
SR_df[SR_df == "."] <- ","
SR_df[SR_df == "BIO1_KH119_A"] <- 1
SR_df[SR_df == "BIO2_KH120_B"] <- 2
SR_df[SR_df == "BIO3_KH121_A"] <- 3
SR_df[SR_df == "BIO4_KH122_B"] <- 4
SR_df[SR_df == "BIO5_KH123_A"] <- 5
SR_df[SR_df == "BIO6_KH124_B"] <- 6

SR_table <- t(SR_df)
colnames(SR_table) <- rownames(SR_df)
SR_table <- SR_table[-1,]

# Decide which rows to keep
to_keep <- c('Mean Reads per Spot', "Median Genes per Spot", "Median UMI Counts per Spot")
SR_table <- SR_table[(row.names(SR_table) %in% to_keep), ]
class(SR_table) <- "numeric"
SR_table <- round(SR_table)

# Create table
SR_table <- tibble::rownames_to_column(as.data.frame(SR_table), "Measurement")
gt(SR_table) |> tab_spanner(label = "Visium sample number", columns = c(2:7)) |>
  tab_options(table.font.names = "Cambria", table.font.size = 12) |>
  gt::gtsave(filename = "Table2_SR_output.png")


# Table 6 ----
setwd("E:/MSc_EMC_project/Main_project/02_BatchCorr_Niches_outs/filtereddata/volcano_plots/")
genes_of_interest <- c("ACTC1", "XIRP2", "CXCL9", "MYH2", "TNNT3", "KRT1")

# Load TCMR vs IF/TA
TCMRvsIFTA_xlsx <- read_excel("01_TCMRvsIFTA.xlsx")
# Add diagnosis
TCMRvsIFTA_xlsx$diag <- "TCMR vs IF/TA"
# Obtain rows with genes of interest
TCMRvsIFTA_xlsx <- filter(TCMRvsIFTA_xlsx, gene %in% genes_of_interest)

# Load AMR vs IF/TA
AMRvsIFTA_xlsx <- read_excel("02_AMRvsIFTA.xlsx")
# Add diagnosis
AMRvsIFTA_xlsx$diag <- "AMR vs IF/TA"
# Obtain rows with genes of interest
AMRvsIFTA_xlsx <- filter(AMRvsIFTA_xlsx, gene %in% genes_of_interest)

# Load TCMR vs AMR
TCMRvsAMR_xlsx <- read_excel("03_TCMRvsAMR.xlsx")
# Add diagnosis
TCMRvsAMR_xlsx$diag <- "TCMR vs AMR"
# Obtain rows with genes of interest
TCMRvsAMR_xlsx <- filter(TCMRvsAMR_xlsx, gene %in% genes_of_interest)

all_genes_of_interest <- rbind(TCMRvsIFTA_xlsx, AMRvsIFTA_xlsx, TCMRvsAMR_xlsx)
all_genes_of_interest <- all_genes_of_interest[,2:8]
all_genes_of_interest <- all_genes_of_interest[,-2:-4]
all_genes_of_interest <- all_genes_of_interest[,-3]
colnames(all_genes_of_interest) <- c("Log2FC", "Gene", "Diagnosis")

all_genes_of_interest1 <- pivot_wider(all_genes_of_interest, names_from = 'Diagnosis', values_from = 'Log2FC')
all_genes_of_interest1 <- all_genes_of_interest1[c(2,3,1,5,4,6),]

gt(all_genes_of_interest1) |>
  tab_options(table.font.names = "Cambria", table.font.size = 12) |>
  fmt_number(decimals = , columns = c(colnames(all_genes_of_interest1[2:4]))) |>
  tab_spanner(label = "Diagnosis comparison", columns = c(2:4)) |>
  tab_style(style = list(
    cell_fill(color = "#e8c5bd"),
    cell_text(style = "italic")),
    locations = cells_body(
      columns = c(`TCMR vs IF/TA`),
      rows = `TCMR vs IF/TA` > 0)) |>
  
  tab_style(style = list(
      cell_fill(color = "#e8c5bd"),
      cell_text(style = "italic")),
    locations = cells_body(
      columns = c(`AMR vs IF/TA`),
      rows = `AMR vs IF/TA` > 0)) |>

  tab_style(style = list(
    cell_fill(color = "#e8c5bd"),
    cell_text(style = "italic")),
    locations = cells_body(
      columns = c(`TCMR vs AMR`),
      rows = `TCMR vs AMR` > 0)) |>
  
  tab_style(style = list(
    cell_fill(color = "#bddee8"),
    cell_text(style = "italic")),
    locations = cells_body(
      columns = c(`TCMR vs IF/TA`),
      rows = `TCMR vs IF/TA` < 0)) |>
  
  tab_style(style = list(
    cell_fill(color = "#bddee8"),
    cell_text(style = "italic")),
    locations = cells_body(
      columns = c(`AMR vs IF/TA`),
      rows = `AMR vs IF/TA` < 0)) |>
  
  tab_style(style = list(
    cell_fill(color = "#bddee8"),
    cell_text(style = "italic")),
    locations = cells_body(
      columns = c(`TCMR vs AMR`),
      rows = `TCMR vs AMR` < 0)) |>
  tab_options(data_row.padding = pct(1.5)) |>
  gt::gtsave(filename = "Table6_DGEA_genesofinterest.png", path = "E:/MSc_EMC_project/Main_project/02_BatchCorr_Niches_outs/filtereddata/")
