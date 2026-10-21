# TyMD-Schistosoma-Infection-Metabolomics

> **Companion repository for:**  
> Kameni Poumeni M. et al. "Asymptomatic *Schistosoma mansoni* infection robustly associates with altered tyrosine metabolism in the plasma of school-aged children."  
> *Nature Communications* (under revision, 2026)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R version](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)
[![MetaboLights](https://img.shields.io/badge/MetaboLights-REQ20260416218755-green)](https://www.ebi.ac.uk/metabolights/)
[![Metabolomics Workbench](https://img.shields.io/badge/Metabolomics%20Workbench-Kameni237__20260425__110834-orange)](https://www.metabolomicsworkbench.org/)

---

## Overview

This repository contains all R analysis scripts, input data tables, and output figures required to reproduce the metabolomics analyses described in the above manuscript. The study identified a **Tyrosine Metabolism-Derived (TyMD) plasma signature** that robustly distinguishes *S. mansoni*-infected school-aged children from uninfected controls across three independent cohorts, using untargeted LC-MS/MS plasma metabolomics.

A companion repository for the hepatic morbidity paper (Paper 2 — ArgPro/FPAS signature) is available at:  
**[ArgPro-Schistosoma-Fibrosis-Metabolomics](https://github.com/DrKameni/ArgPro-Schistosoma-Fibrosis-Metabolomics)**

---

## Study Design

### Populations

| Cohort | Design | KK+ Infected | KK− Control | Platform | Notes |
|--------|--------|:---:|:---:|----------|-------|
| **Discovery** | Cross-sectional | 38 | 31 | ZIC-HILIC Orbitrap (Glasgow Polyomics) | 7/10 KK− controls CAA-positive (sub-patent); included as controls per pre-specified protocol |
| **Validation** | Cross-sectional | 32 | 28 | ZIC-HILIC Orbitrap (Glasgow Polyomics) | Independent school, same region |
| **Polyparasitism** | Cross-sectional | 29 | 26 | ZIC-HILIC Orbitrap (Glasgow Polyomics) | Co-infected with soil-transmitted helminths |

**Total N:** 184 participants (99 KK+, 85 KK−)  
**Age range:** 5–17 years (school-aged children)  
**Setting:** Cameroon (sub-Saharan Africa)  
**Infection criterion:** Kato-Katz (KK) stool microscopy ≥ 1 egg/slide

> ⚠️ **Critical note on Discovery controls:** Seven of the 31 KK− controls in the Discovery cohort were CAA (Circulating Anodic Antigen) positive, indicating sub-patent infection. These were retained as controls per the pre-specified study protocol. Sensitivity analyses excluding CAA-positive individuals are reported in Supplementary Table S8.

### Three-Phase Analytical Design

| Phase | Cohort | Metabolite selection | ML method | Validation approach |
|-------|--------|---------------------|-----------|---------------------|
| **Phase 1: Discovery** | Discovery (n=69) | Volcano (FC≥1.2, p<0.1) | RF LOOCV | Internal LOOCV only |
| **Phase 2: Validation** | Validation (n=60) | Volcano (FC≥1.3, p<0.05) | RF LOOCV | Independent cohort |
| **Phase 3: Polyparasitism** | Polyparasitism (n=55) | Volcano (FC≥1.5, p<0.05) | RF LOOCV | Co-infection generalizability |

**Volcano FC thresholds** differ across phases to account for increasing biological noise with parasitic co-infection burden.

---

## Key Findings

### Machine Learning Performance (RF LOOCV)

| Cohort | AUC (95% CI) | Sensitivity | Specificity |
|--------|-------------|:-----------:|:-----------:|
| Discovery | 0.885 (0.715–1.000) | 0.868 | 0.871 |
| Validation | 0.904 (0.812–0.996) | 0.875 | 0.893 |
| Polyparasitism | 0.971 (0.922–1.000) | 0.966 | 0.923 |

### Meta-Analysis (Random Effects, k=3 batches)

| Model | Pooled AUC | 95% CI | I² | Q-test p |
|-------|-----------|--------|:--:|:--------:|
| DerSimonian-Laird | **0.911** | 0.823–0.960 | 0.0% | 0.847 |

### Tyrosine Pathway Activity Score (TPAS)

The TPAS is a direction-corrected composite z-score computed within each batch:

- **Wilcoxon p-values:** Discovery p=0.008, Validation p=0.002, Polyparasitism p<0.001
- **Combined-cohort LOOCV AUC:** 0.847 (95% CI: 0.786–0.908)
- **TPAS correlates with EPG** (Spearman ρ=0.41, p<0.001) — higher egg burden = higher TPAS

**TyMD signature features** (shared across ≥2 cohorts):

| Feature | Cohorts | Direction (KK+ vs KK−) |
|---------|---------|------------------------|
| Tyrosine | Disc, Val | ↑ Up |
| Tyramine | All 3 | ↑ Up |
| Hydroxyphenylacetylglycine | Disc, Poly | ↑ Up |
| 4-Hydroxycyclohexylcarboxylic acid | Val, Poly | ↑ Up |
| (4-Hydroxy-3-methoxyphenyl)ethanol | Poly | ↑ Up |
| 3-Methoxytyramine | Disc, Val | ↑ Up |

> ⚠️ **Leakage disclaimer:** RF LOOCV AUC values reported above are **optimistic upper bounds** because volcano-based feature selection was performed on the **full batch** before LOOCV, constituting feature-selection leakage. True generalisation performance is better estimated by the **meta-analysis pooled AUC** or the **combined-cohort LOOCV AUC**. This is acknowledged in the manuscript and reviewer response (MDF comment #59).

---

## Repository Structure

```
TyMD-Schistosoma-Infection-Metabolomics/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── 00_Master_Reproducibility.R        # Orchestrates all steps in sequence
│
├── scripts/
│   ├── 01_Boxplots_Tyrosine.R         # Tyrosine pathway DEM boxplots (3 batches)
│   ├── 02_ML_RandomForest_LOOCV.R     # RF LOOCV + SHAP importance (all batches)
│   ├── 03_TPAS_Computation.R          # TPAS score per batch (direction-corrected z-score)
│   ├── 04_MetaAnalysis_AUC.R          # Random-effects meta-analysis (DL, k=3)
│   ├── 05_TPAS_CombinedCohort.R       # Combined-cohort LOOCV on TPAS
│   ├── 06_Clinical_Correlations.R     # TPAS vs EPG, spleen, liver (Spearman + BH-FDR)
│   ├── 07_ChemClass_Waffle.R          # Chemical class waffle charts per batch
│   ├── 08_DiscPower_Ranking.R         # Discriminatory power of metabolite categories
│   └── 09_Sensitivity_Analysis.R      # CAA-positive exclusion sensitivity analysis
│
├── data/
│   ├── Discovery_ML_input.csv         # n=69; feature matrix + Group label
│   ├── Validation_ML_input.csv        # n=60; feature matrix + Group label
│   ├── Polyparasitism_ML_input.csv    # n=55; feature matrix + Group label
│   ├── DEM_annotation_reference.csv   # Metabolite annotations (HMDB IDs, KEGG IDs)
│   ├── Clinical_master.csv            # Clinical variables (EPG, spleen score, age, sex)
│   └── TPAS_scores_allbatches.csv     # Output of Step 03 (used by Steps 04–06)
│
├── figures/                           # All manuscript figures (PDF + PNG)
│   ├── Fig1_StudyDesign.pdf
│   ├── Fig2_VolcanoPlots.pdf
│   ├── Fig3_Boxplots_TyrPathway.pdf
│   ├── Fig4_RFLOOCV_ROC.pdf
│   ├── Fig5_SHAP_Importance.pdf
│   ├── Fig6_TPAS_scores.pdf
│   ├── Fig7_MetaAnalysis_Forest.pdf
│   └── FigS_ChemClass_Waffle.pdf
│
└── tables/
    ├── Table1_ParticipantCharacteristics.csv
    ├── Table2_RF_Performance.csv
    ├── Table3_SHAP_Top15.csv
    └── TableS_TyrFeatures_RF_SHAP_Ranks.csv
```

---

## Installation

### R Dependencies

```r
# Install all required packages
required_packages <- c(
  # Core analysis
  "tidyverse", "readxl", "writexl",
  # Machine learning
  "randomForest", "iml", "pROC",
  # Statistics
  "meta", "coin",
  # Visualization
  "ggplot2", "ggbeeswarm", "waffle", "ggrepel", "patchwork",
  # Tables
  "gt", "flextable"
)

# Install missing packages
missing_pkgs <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

# waffle chart package (from GitHub if not on CRAN)
if (!"waffle" %in% installed.packages()[,"Package"]) {
  remotes::install_github("hrbrmstr/waffle")
}
```

### R Version

Analyses were performed with **R version 4.3.x**. Key package versions:

| Package | Version tested |
|---------|---------------|
| randomForest | 4.7-1.1 |
| iml | 0.11.1 |
| pROC | 1.18.5 |
| meta | 7.0-0 |
| tidyverse | 2.0.0 |

---

## Usage

### Option 1 — Run all steps sequentially (recommended)

```r
# Set your working directory to the repository root
setwd("/path/to/TyMD-Schistosoma-Infection-Metabolomics")

# Source the master script (runs all 9 steps in order)
source("00_Master_Reproducibility.R")
```

### Option 2 — Run individual steps

```r
# Example: run only the machine learning step
source("scripts/02_ML_RandomForest_LOOCV.R")
```

> Steps 04–06 depend on `data/TPAS_scores_allbatches.csv` produced by Step 03.  
> Run Step 03 before Steps 04–06.

### Reproducibility

- **Random seed:** `set.seed(42)` is applied at the start of every script that uses random sampling (RF LOOCV, permutation tests)
- **Missing value imputation:** Half-minimum per column (left-censored LC-MS/MS convention); applied independently within each batch
- **Volcano thresholds:** Discovery FC≥1.2, p<0.1; Validation FC≥1.3, p<0.05; Polyparasitism FC≥1.5, p<0.05
- **FDR correction:** Benjamini-Hochberg (BH) throughout
- **Expected runtime:** ~25–40 min on a standard laptop (LOOCV with SHAP is the bottleneck)

---

## Statistical Methods Summary

| Step | Script | Method | Key output |
|------|--------|--------|-----------|
| 1 | `01_Boxplots_Tyrosine.R` | Wilcoxon rank-sum + BH-FDR | Boxplots per Tyr pathway metabolite |
| 2 | `02_ML_RandomForest_LOOCV.R` | RF LOOCV (ntree=1000) + SHAP | AUC, ROC curves, SHAP importances |
| 3 | `03_TPAS_Computation.R` | Direction-corrected composite z-score | `TPAS_scores_allbatches.csv` |
| 4 | `04_MetaAnalysis_AUC.R` | DerSimonian-Laird random effects | Pooled AUC, I², forest plot |
| 5 | `05_TPAS_CombinedCohort.R` | Logistic LOOCV on pooled TPAS | AUC, sensitivity, specificity |
| 6 | `06_Clinical_Correlations.R` | Spearman ρ + BH-FDR | TPAS vs EPG/spleen/liver heatmap |
| 7 | `07_ChemClass_Waffle.R` | Chemical class annotation | Waffle charts per batch |
| 8 | `08_DiscPower_Ranking.R` | AUC per metabolite category | Bubble plot + ranking bar |
| 9 | `09_Sensitivity_Analysis.R` | Exclusion of CAA+ controls | Sensitivity AUC (Discovery) |

---

## Data Availability

Raw metabolomics data have been deposited in two public repositories:

| Repository | Accession / ID | URL |
|-----------|---------------|-----|
| **MetaboLights** (EBI) | REQ20260416218755 | https://www.ebi.ac.uk/metabolights/ |
| **Metabolomics Workbench** | Kameni237_20260425_110834 | https://www.metabolomicsworkbench.org/ |

The input CSV files in `data/` are derived from the PiMP pipeline (Glasgow Polyomics) and correspond directly to the deposited datasets.

---

## Citing This Work

If you use these scripts or data, please cite:

> Kameni Poumeni, M. et al. Asymptomatic *Schistosoma mansoni* infection robustly associates with altered tyrosine metabolism in the plasma of school-aged children. *Nature Communications* (2026). https://github.com/DrKameni/TyMD-Schistosoma-Infection-Metabolomics

---

## License

This repository is released under the [MIT License](LICENSE).

---

## Contact

**Dr. Mireille Kameni Poumeni**  
Email: my.kameni@gmail.com  
GitHub: [@DrKameni](https://github.com/DrKameni)

---

*Last updated: October 2026*
