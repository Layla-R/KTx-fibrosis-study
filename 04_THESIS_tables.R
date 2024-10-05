# Libraries ----
library(dplyr) 
library(gt)
library(gtsummary)
library(readxl)

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

# Table Figure 3 - dotplot ----
cluster_no <- c(as.character(0:10))
annotation <- c("loop of Henle",
                "Proximal tubules",
                "Distal tubules",
                "???",
                "Transmembrane transporter, ???",
                "Mitochondrial",
                "Ceccular response/renal homeostasis",
                "B-cell/antigen receptor pathway",
                "Glomeruli",
                "Fibrotic",
                "Muscle contraction/myofibril, ??")

annotation_table <- data.frame(cluster_no, annotation)
colnames(annotation_table) <- c("Cluster", "Annotation")

gt(annotation_table) |> 
  tab_options(table.font.names = "Cambria", table.font.size = 12)
