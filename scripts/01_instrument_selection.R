# ----------------------------------------
# Step 1: Instrument Selection
# Author: Neha Manoti
#
# Description:
# Select genome-wide significant SNPs from proteomics GWAS
# and filter based on instrument strength (F-statistic > 10)
#
# Output:
# - formatted_exposures.rds (processed SNP instruments per protein)
# ----------------------------------------


# Load libraries
library(data.table)
library(stringr)
library(dplyr)

##Define project-relative paths (ensures reproducibility across systems)
base_dir <- "."

exposure_location <- file.path(base_dir, "data", "exposure")

output_dir <- file.path(base_dir, "data", "processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List exposure files
exposure_files <- list.files(
  exposure_location,
  pattern = "\\.csv$",
  full.names = TRUE
)

# SNP formatting function
format_snpids <- function(df, id_column){
  
  df[[id_column]] <- as.character(df[[id_column]])
  
  pat_full <- "^X?([^\\.]+)\\.([0-9]+)_([^_]+)_([^_]+)_([^_]+)$"
  m <- stringr::str_match(df[[id_column]], pat_full)
  
  df$snp_id <- NA_character_
  
  matched <- which(!is.na(m[,1]))
  
  if(length(matched) > 0){
    chr <- m[matched, 2]
    pos <- m[matched, 3]
    ref <- m[matched, 4]
    alt <- m[matched, 5]
    
    df$snp_id[matched] <- paste0(chr, ":", pos, "_", ref, "_", alt)
  }
  
  return(df)
}

# Process all exposure files
formatted_exposures <- list()

for(f in exposure_files){
  
  exposure_data <- fread(f)
  
  exposure_data <- format_snpids(exposure_data, "rsid")
  
  exposure_data <- exposure_data[, .(
    SNP = snp_id,
    beta = beta,
    se = SE,
    effect_allele = EA,
    other_allele = OA,
    eaf = EAF,
    pval = p,
    chr = chr,
    pos = pos
  )]
  
  # Remove missing values
  exposure_data <- exposure_data[complete.cases(exposure_data), ]
  
  # Calculate F-statistic
  exposure_data$f_stat <- (exposure_data$beta^2) / (exposure_data$se^2)
  
  # Keep strong instruments
  exposure_data <- exposure_data[exposure_data$f_stat > 10, ]
  
  protein_name <- tools::file_path_sans_ext(basename(f))
  
  formatted_exposures[[protein_name]] <- exposure_data
}

# Save processed exposures
saveRDS(
  formatted_exposures,
  file = file.path(output_dir, "formatted_exposures.rds")
)

# Summary stats
cat("Number of proteins:", length(formatted_exposures), "\n")
