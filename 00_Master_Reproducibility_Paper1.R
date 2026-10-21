##############################################################################
#  00_Master_Reproducibility_Paper1.R
#
#  Repository : TyMD-Schistosoma-Infection-Metabolomics
#  Paper      : "Asymptomatic Schistosoma mansoni infection robustly
#               associates with altered tyrosine metabolism in the plasma
#               of school-aged children."
#  Authors    : Kameni Poumeni M. et al. | Nature Communications (2026)
#
#  Purpose    : Run all 9 analysis steps in order, check dependencies,
#               verify expected outputs, and print a reproducibility report.
#
#  Usage      : setwd("/path/to/TyMD-Schistosoma-Infection-Metabolomics")
#               source("00_Master_Reproducibility_Paper1.R")
#
#  Runtime    : ~25-40 min on a standard laptop (RF LOOCV bottleneck)
#  R version  : >= 4.2.0 recommended
##############################################################################

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║   TyMD-Schistosoma-Infection-Metabolomics  —  Reproducibility Run   ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")
cat("R version:", R.version.string, "\n")
cat("Date     :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Global seed (applied here AND at the start of each script)
# ─────────────────────────────────────────────────────────────────────────────
set.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Package check
# ─────────────────────────────────────────────────────────────────────────────
required_packages <- c(
  "tidyverse", "readxl", "writexl",
  "randomForest", "iml", "pROC",
  "meta", "coin",
  "ggplot2", "ggbeeswarm", "waffle", "ggrepel", "patchwork",
  "gt", "flextable"
)

cat("─── Checking required packages ──────────────────────────────────────────\n")
missing_pkgs <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse=", "), "\n")
  install.packages(missing_pkgs, repos="https://cloud.r-project.org")
  if ("waffle" %in% missing_pkgs) {
    if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
    remotes::install_github("hrbrmstr/waffle")
  }
} else {
  cat("✓ All required packages are installed.\n\n")
}

suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only=TRUE))
})

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Input data file check
# ─────────────────────────────────────────────────────────────────────────────
cat("─── Checking input data files ───────────────────────────────────────────\n")

required_files <- c(
  "data/Discovery_ML_input.csv",
  "data/Validation_ML_input.csv",
  "data/Polyparasitism_ML_input.csv",
  "data/DEM_annotation_reference.csv",
  "data/Clinical_master.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  cat("⚠ Missing input files:\n")
  for (f in missing_files) cat("  -", f, "\n")
  cat("\nDownload raw data from:\n")
  cat("  MetaboLights : https://www.ebi.ac.uk/metabolights/ (REQ20260416218755)\n")
  cat("  Metabolomics Workbench : Kameni237_20260425_110834\n\n")
  stop("Please ensure all input data files are present before proceeding.")
} else {
  cat("✓ All input data files found.\n\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Helper — run a step and time it
# ─────────────────────────────────────────────────────────────────────────────
run_step <- function(step_num, script_name, description) {
  cat(sprintf("─── Step %02d: %s ─────────────────────────────────────────\n",
              step_num, description))
  t0 <- proc.time()
  tryCatch({
    source(file.path("scripts", script_name), local=FALSE)
    elapsed <- (proc.time() - t0)["elapsed"]
    cat(sprintf("✓ Step %02d completed in %.1f sec\n\n", step_num, elapsed))
    return(list(ok=TRUE, time=elapsed))
  }, error=function(e) {
    cat(sprintf("✗ Step %02d FAILED: %s\n\n", step_num, conditionMessage(e)))
    return(list(ok=FALSE, time=NA, error=conditionMessage(e)))
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Run all steps
# ─────────────────────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  RUNNING ANALYSIS PIPELINE\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

results <- list()

results[[1]]  <- run_step(1,  "01_Boxplots_Tyrosine.R",
                           "Tyrosine pathway DEM boxplots")

results[[2]]  <- run_step(2,  "02_ML_RandomForest_LOOCV.R",
                           "RF LOOCV + SHAP (all batches)")

results[[3]]  <- run_step(3,  "03_TPAS_Computation.R",
                           "TPAS score computation (direction-corrected z-score)")

# Steps 4-6 require TPAS_scores_allbatches.csv from Step 3
if (!file.exists("data/TPAS_scores_allbatches.csv")) {
  cat("⚠ data/TPAS_scores_allbatches.csv not found — Steps 4-6 cannot run.\n")
  cat("  Ensure Step 3 completed successfully.\n\n")
  results[[4]] <- list(ok=FALSE, time=NA, error="Dependency missing: TPAS_scores_allbatches.csv")
  results[[5]] <- list(ok=FALSE, time=NA, error="Dependency missing: TPAS_scores_allbatches.csv")
  results[[6]] <- list(ok=FALSE, time=NA, error="Dependency missing: TPAS_scores_allbatches.csv")
} else {
  results[[4]]  <- run_step(4,  "04_MetaAnalysis_AUC.R",
                             "Random-effects meta-analysis (DL, k=3)")

  results[[5]]  <- run_step(5,  "05_TPAS_CombinedCohort.R",
                             "Combined-cohort LOOCV on pooled TPAS")

  results[[6]]  <- run_step(6,  "06_Clinical_Correlations.R",
                             "TPAS vs clinical variables (Spearman + BH-FDR)")
}

results[[7]]  <- run_step(7,  "07_ChemClass_Waffle.R",
                           "Chemical class waffle charts")

results[[8]]  <- run_step(8,  "08_DiscPower_Ranking.R",
                           "Discriminatory power ranking by metabolite category")

results[[9]]  <- run_step(9,  "09_Sensitivity_Analysis.R",
                           "Sensitivity analysis (CAA+ control exclusion)")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Expected key results (for verification)
# ─────────────────────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  EXPECTED KEY RESULTS\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat("RF LOOCV AUC (feature-selection leakage caveat applies):\n")
cat("  Discovery     : 0.885 (95% CI: 0.715-1.000)\n")
cat("  Validation    : 0.904 (95% CI: 0.812-0.996)\n")
cat("  Polyparasitism: 0.971 (95% CI: 0.922-1.000)\n\n")

cat("Meta-analysis (DerSimonian-Laird, k=3 batches):\n")
cat("  Pooled AUC    : 0.911 (95% CI: 0.823-0.960)\n")
cat("  I²            : 0.0%\n")
cat("  Q-test p-value: 0.847\n\n")

cat("TPAS Wilcoxon test (KK+ vs KK-):\n")
cat("  Discovery     : p = 0.008\n")
cat("  Validation    : p = 0.002\n")
cat("  Polyparasitism: p < 0.001\n\n")

cat("Combined-cohort LOOCV (logistic model on pooled TPAS):\n")
cat("  AUC           : 0.847 (95% CI: 0.786-0.908)\n\n")

cat("TPAS vs EPG (Spearman correlation, all batches pooled):\n")
cat("  rho = 0.41, p < 0.001\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Pipeline summary
# ─────────────────────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  PIPELINE SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

step_names <- c(
  "Tyrosine boxplots",
  "RF LOOCV + SHAP",
  "TPAS computation",
  "Meta-analysis",
  "Combined LOOCV",
  "Clinical correlations",
  "Chemical class waffle",
  "Discriminatory power",
  "Sensitivity analysis"
)

n_ok   <- sum(sapply(results, function(x) isTRUE(x$ok)))
n_fail <- length(results) - n_ok
total_time <- sum(sapply(results, function(x) ifelse(is.na(x$time), 0, x$time)))

for (i in seq_along(results)) {
  status <- if (isTRUE(results[[i]]$ok)) "✓" else "✗"
  t_str  <- if (!is.na(results[[i]]$time)) sprintf("%.1f sec", results[[i]]$time) else "—"
  cat(sprintf("  %s  Step %02d  %-30s  %s\n",
              status, i, step_names[i], t_str))
  if (!isTRUE(results[[i]]$ok) && !is.null(results[[i]]$error)) {
    cat(sprintf("           ERROR: %s\n", results[[i]]$error))
  }
}

cat(sprintf("\nSteps completed: %d / %d   |   Total runtime: %.1f sec (%.1f min)\n",
            n_ok, length(results), total_time, total_time/60))

if (n_fail == 0) {
  cat("\n✓ All steps completed successfully. Figures in figures/ — Tables in tables/\n")
  cat("  Cite as: Kameni Poumeni M. et al. Nature Communications (2026)\n")
} else {
  cat(sprintf("\n⚠ %d step(s) failed. Check error messages above.\n", n_fail))
}

cat("\n")
cat("──────────────────────────────────────────────────────────────────────\n")
cat("  Data deposited at:\n")
cat("  MetaboLights          : REQ20260416218755\n")
cat("  Metabolomics Workbench: Kameni237_20260425_110834\n")
cat("──────────────────────────────────────────────────────────────────────\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7.  Session info (for reproducibility logging)
# ─────────────────────────────────────────────────────────────────────────────
cat("─── Session information ─────────────────────────────────────────────────\n")
print(sessionInfo())
