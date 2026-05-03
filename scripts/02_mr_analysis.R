# ----------------------------------------
# Step 2: Outcome Formatting
# Author: Neha Manoti
#
# Description:
# This script formats outcome GWAS summary statistics (e.g. BAVS)
# for compatibility with the TwoSampleMR package.
#
# Key steps:
# - Reads outcome GWAS file
# - Standardizes column names
# - Converts to MR-compatible format
#
# Output:
# - outcome_formatted.rds
# ----------------------------------------

# Load libraries
library(data.table)
library(TwoSampleMR)

# Define paths
base_dir <- "."

outcome_file <- file.path(base_dir, "data", "outcome", "outcome_data.tsv.gz")

# NOTE:
# Replace with your GWAS outcome file.
# Data not included due to access restrictions.

output_dir <- file.path(base_dir, "data", "processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load outcome data
outcome_data <- fread(outcome_file)

# Format for TwoSampleMR
outcome_formatted <- format_data(
  outcome_data,
  type = "outcome",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  chr_col = "chromosome",
  pos_col = "base_pair_location",
  samplesize_col = "n"
)

# Save formatted outcome
saveRDS(
  outcome_formatted,
  file = file.path(output_dir, "outcome_formatted.rds")
)

# Quick check
cat("Number of SNPs in outcome:", nrow(outcome_formatted), "\n")
