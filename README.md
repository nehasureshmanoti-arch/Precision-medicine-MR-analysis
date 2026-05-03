🧬 Mendelian Randomization Pipeline for Target Identification

Author: Neha Manoti
Affiliation: Uppsala University
Project: MSc Precision Medicine Thesis

📌 Overview

This repository contains a reproducible pipeline for performing two-sample Mendelian Randomization (MR) to investigate causal relationships between circulating proteins and disease outcomes.

The workflow integrates large-scale proteomics and GWAS datasets (UK Biobank PPP) to identify potential therapeutic targets and prioritise candidates for downstream validation.

🎯 Objective

To identify proteins with causal effects on disease risk using genetic instruments and robust MR methodology, supporting target discovery in precision medicine.

## 🔁 Workflow

- **Instrument Selection**  
  - Selection of genome-wide significant SNPs  
  - Filtering based on instrument strength (F-statistic > 10)

- **Outcome Formatting**  
  - Standardisation of GWAS summary statistics  
  - Conversion to MR-compatible format  

- **Harmonisation**  
  - Alignment of exposure and outcome datasets  
  - Removal of ambiguous or mismatched SNPs  

- **Mendelian Randomization Analysis**  
  - Inverse Variance Weighted (IVW)  
  - MR-Egger regression  
  - Weighted median method  
  - Sensitivity analyses (heterogeneity, pleiotropy)  

- **Visualisation**  
  - Forest plots for effect estimates  
  - Bubble lattice plots across MR methods

## 📁 Repository Structure

```text
scripts/
  01_instrument_selection.R
  02_format_outcome.R
  03_harmonisation.R
  04_mr_analysis.R
  05_visualization.R

data/
  exposure/
  outcome/
  processed/

results/
  plots/

docs/
  methods.md

## ▶️ How to Run

1. Place exposure and outcome GWAS files in the `data/` directory  
2. Run scripts in order:

- scripts/01_instrument_selection.R  
- scripts/02_format_outcome.R  
- scripts/03_harmonisation.R  
- scripts/04_mr_analysis.R  
- scripts/05_visualization.R  

3. Results will be saved in the `results/` folder
  
## 📊 Key Outputs

- MR results table (`MR_results_table.csv`)  
- Heterogeneity and pleiotropy results  
- Filtered robust associations (`robust_MR_hits.csv`)  
- Publication-ready plots (forest plot, bubble lattice)  

---

## 🔒 Data Availability

Due to data access restrictions (e.g. UK Biobank), raw datasets are not included in this repository.

Users should place their input files in the appropriate `data/` directories and update file paths if necessary.

---

## 🧪 Tools & Technologies

- R (TwoSampleMR, data.table, dplyr, ggplot2)  
- GWAS summary statistics  
- Olink proteomics (UK Biobank PPP)  

---

## 📌 Notes

This repository is designed to ensure **transparency, reproducibility, and modularity** of Mendelian Randomization analyses performed during my MSc thesis.
