# Libraries ----
library(dplyr) 
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

