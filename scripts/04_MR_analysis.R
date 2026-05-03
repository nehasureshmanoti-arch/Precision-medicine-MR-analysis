# ----------------------------------------
# Step 4: Mendelian Randomization Analysis
# Author: Neha Manoti
#
# Description:
# Performs MR analysis using multiple methods and evaluates
# robustness through sensitivity tests.
#
# Key steps:
# - IVW, MR-Egger, Weighted Median
# - Heterogeneity testing
# - Pleiotropy assessment
# - Basic filtering of robust associations
#
# Output:
# - MR_results_table.csv
# - heterogeneity_results.csv
# - pleiotropy_results.csv
# - robust_MR_hits.csv
# ----------------------------------------

# Load libraries
library(TwoSampleMR)
library(data.table)
library(dplyr)

# Define paths
base_dir <- "."

harmonised_file <- file.path(base_dir, "data", "processed", "harmonised_data.rds")
output_dir <- file.path(base_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
harmonised_list <- readRDS(harmonised_file)

# -----------------------------
# MR ANALYSIS
# -----------------------------
mr_results <- list()

for (prot in names(harmonised_list)) {
  
  df <- harmonised_list[[prot]]
  
  if (nrow(df) >= 1) {
    
    res <- mr(
      df,
      method_list = c(
        "mr_ivw",
        "mr_egger_regression",
        "mr_weighted_median"
      )
    )
    
    if (nrow(res) > 0) {
      res$protein <- prot
      mr_results[[prot]] <- res
    }
  }
}

mr_results_df <- do.call(rbind, mr_results)

write.csv(
  mr_results_df,
  file.path(output_dir, "MR_results_table.csv"),
  row.names = FALSE
)

# -----------------------------
# HETEROGENEITY
# -----------------------------
heterogeneity_results <- list()

for (prot in names(harmonised_list)) {
  
  dat <- harmonised_list[[prot]]
  
  if (nrow(dat) >= 3) {
    
    q <- mr_heterogeneity(dat)
    
    if (nrow(q) > 0) {
      q$protein <- prot
      heterogeneity_results[[prot]] <- q
    }
  }
}

heterogeneity_df <- do.call(rbind, heterogeneity_results)

write.csv(
  heterogeneity_df,
  file.path(output_dir, "heterogeneity_results.csv"),
  row.names = FALSE
)

# -----------------------------
# PLEIOTROPY
# -----------------------------
pleiotropy_results <- list()

for (prot in names(harmonised_list)) {
  
  dat <- harmonised_list[[prot]]
  
  if (nrow(dat) >= 3) {
    
    p <- mr_pleiotropy_test(dat)
    
    if (nrow(p) > 0) {
      p$protein <- prot
      pleiotropy_results[[prot]] <- p
    }
  }
}

pleiotropy_df <- do.call(rbind, pleiotropy_results)

write.csv(
  pleiotropy_df,
  file.path(output_dir, "pleiotropy_results.csv"),
  row.names = FALSE
)

# -----------------------------
# FILTER ROBUST HITS (SIMPLE VERSION)
# -----------------------------
mr_ivw <- mr_results_df %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(pval)

mr_ivw_sig <- mr_ivw %>%
  filter(pval < 0.05, nsnp >= 3)

# Remove heterogeneity hits
bad_het <- heterogeneity_df %>%
  filter(Q_pval < 0.05, method == "Inverse variance weighted") %>%
  pull(protein)

mr_ivw_clean <- mr_ivw_sig %>%
  filter(!protein %in% bad_het)

write.csv(
  mr_ivw_clean,
  file.path(output_dir, "robust_MR_hits.csv"),
  row.names = FALSE
)

# Summary
cat("Total MR results:", nrow(mr_results_df), "\n")
cat("Robust hits:", nrow(mr_ivw_clean), "\n")
