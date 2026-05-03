# ----------------------------------------
# Step 3: Harmonisation
# Author: Neha Manoti
#
# Description:
# Aligns exposure and outcome datasets to ensure consistent
# effect allele orientation prior to Mendelian Randomization.
#
# Key steps:
# - Load processed exposure data
# - Load formatted outcome data
# - Convert exposure to MR format
# - Harmonise datasets using TwoSampleMR
#
# Output:
# - harmonised_data.rds
# ----------------------------------------

# Load libraries
library(TwoSampleMR)
library(data.table)

# Define paths
base_dir <- "."

exposure_file <- file.path(base_dir, "data", "processed", "formatted_exposures.rds")
outcome_file  <- file.path(base_dir, "data", "processed", "outcome_formatted.rds")

output_dir <- file.path(base_dir, "data", "processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
formatted_exposures <- readRDS(exposure_file)
outcome_formatted   <- readRDS(outcome_file)

# Initialise list
harmonised_list <- list()

# Loop through each protein
for (prot in names(formatted_exposures)) {
  
  exp_raw <- formatted_exposures[[prot]]
  
  # Skip empty datasets
  if (nrow(exp_raw) == 0) next
  
  exp_df <- as.data.frame(exp_raw)
  
  # Format exposure for MR
  exp_fmt <- format_data(
    exp_df,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
  )
  
  exp_fmt$exposure <- prot
  
  # Harmonise
  harm <- harmonise_data(
    exposure_dat = exp_fmt,
    outcome_dat  = outcome_formatted
  )
  
  # Store only if SNPs remain
  if (nrow(harm) > 0) {
    harmonised_list[[prot]] <- harm
  }
}

# Save harmonised data
saveRDS(
  harmonised_list,
  file = file.path(output_dir, "harmonised_data.rds")
)

# Summary
cat("Number of proteins after harmonisation:", length(harmonised_list), "\n")
