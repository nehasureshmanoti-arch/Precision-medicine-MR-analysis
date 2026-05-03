# ----------------------------------------
# Step 5: Visualization
# Author: Neha Manoti
#
# Description:
# Generates plots to visualise MR results and highlight
# significant protein–disease associations.
#
# Output:
# - MR_forestplot.pdf
# - bubble_lattice.pdf
# ----------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)

# Paths
base_dir <- "."
results_dir <- file.path(base_dir, "results")
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Load MR results
mr_results <- read.csv(file.path(results_dir, "MR_results_table.csv"))

# -----------------------------
# FOREST PLOT (IVW only)
# -----------------------------
mr_ivw <- mr_results %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(pval)

mr_ivw_sig <- mr_ivw %>%
  filter(pval < 0.05, nsnp >= 3)

if (nrow(mr_ivw_sig) > 0) {
  
  mr_ivw_sig <- mr_ivw_sig %>%
    mutate(
      lower = b - 1.96 * se,
      upper = b + 1.96 * se
    )
  
  p_forest <- ggplot(mr_ivw_sig,
                     aes(x = reorder(protein, b), y = b)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_classic() +
    labs(
      x = "Protein",
      y = "MR effect estimate (IVW)",
      title = "Proteins associated with outcome"
    )
  
  ggsave(
    file.path(plots_dir, "MR_forestplot.pdf"),
    p_forest,
    width = 7,
    height = max(4, 0.35 * nrow(mr_ivw_sig)),
    dpi = 300
  )
}

# -----------------------------
# BUBBLE LATTICE (TOP PROTEINS)
# -----------------------------
mr_plot <- mr_results %>%
  filter(method %in% c(
    "Inverse variance weighted",
    "Weighted median",
    "MR Egger"
  ))

# Top proteins from IVW
top_proteins <- mr_plot %>%
  filter(method == "Inverse variance weighted", pval < 0.05, nsnp >= 3) %>%
  arrange(pval) %>%
  slice(1:10) %>%
  pull(protein)

mr_plot <- mr_plot %>%
  filter(protein %in% top_proteins) %>%
  mutate(
    OR = exp(b),
    significance = ifelse(pval < 0.05, "Significant", "Not Significant")
  )

if (nrow(mr_plot) > 0) {
  
  p_bubble <- ggplot(mr_plot, aes(x = method, y = protein)) +
    geom_point(aes(size = OR, color = significance), alpha = 0.8) +
    scale_color_manual(values = c("grey60", "red")) +
    scale_size_continuous(range = c(3, 10)) +
    theme_minimal() +
    labs(
      title = "MR results across methods (top proteins)",
      x = "MR Method",
      y = "Protein",
      size = "Odds Ratio",
      color = "Significance"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(
    file.path(plots_dir, "bubble_lattice.pdf"),
    p_bubble,
    width = 10,
    height = 6
  )
}

cat("Plots generated in:", plots_dir, "\n")
