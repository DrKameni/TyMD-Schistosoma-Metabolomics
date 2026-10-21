##############################################################################
#  09_Sensitivity_Analysis.R
#
#  Sensitivity analysis: exclusion of CAA-positive KK− controls
#  in the Discovery batch.
#
#  Background:
#    7 of the 31 KK− controls in Discovery were CAA (Circulating Anodic
#    Antigen) positive, indicating sub-patent S. mansoni infection.
#    These were retained as controls per the pre-specified study protocol
#    (primary analysis). This script repeats the RF-LOOCV and TPAS
#    analysis after excluding those 7 samples to test robustness.
#
#  Input  : data/Discovery_ML_input.csv
#            data/Clinical_master.csv  (must contain CAA_status column)
#  Output : figures/FigS_Sensitivity_CAA.pdf
#            tables/TableS_Sensitivity_CAA.csv
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(pROC)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ─── Load data ────────────────────────────────────────────────────────────────
disc <- read_csv("data/Discovery_ML_input.csv",  show_col_types = FALSE)
clin <- read_csv("data/Clinical_master.csv",      show_col_types = FALSE)
disc$Group <- factor(disc$Group, levels = c("Control", "Infected"))

cat("Discovery: n =", nrow(disc), "\n")
cat("Groups:\n"); print(table(disc$Group))

# ─── Check CAA status availability ───────────────────────────────────────────
if (!"CAA_status" %in% colnames(clin)) {
  cat("\n⚠ Column 'CAA_status' not found in Clinical_master.csv.\n")
  cat("  Expected values: 'positive' or 'negative'.\n")
  cat("  If CAA data is stored under a different column name, edit this script.\n\n")
  cat("  FALLBACK: using placeholder IDs for the 7 known CAA+ controls.\n")
  cat("  Replace 'caa_ids_placeholder' below with the actual SampleIDs.\n\n")

  # Placeholder — replace with actual SampleIDs of CAA+ controls
  caa_ids_placeholder <- character(0)   # ← insert SampleIDs here if known
  caa_positive_ids <- caa_ids_placeholder
} else {
  caa_df <- left_join(disc %>% select(SampleID), clin, by = "SampleID")
  caa_positive_ids <- caa_df$SampleID[
    caa_df$CAA_status == "positive" & disc$Group == "Control"]
  cat(sprintf("CAA-positive controls identified: n=%d\n", length(caa_positive_ids)))
  cat("IDs:", paste(caa_positive_ids, collapse=", "), "\n\n")
}

# ─── Sensitivity dataset (CAA+ excluded) ─────────────────────────────────────
if (length(caa_positive_ids) == 0) {
  cat("No CAA-positive IDs to exclude — sensitivity = primary analysis.\n")
  disc_sens <- disc
} else {
  disc_sens <- filter(disc, !SampleID %in% caa_positive_ids)
  cat(sprintf("After exclusion: n=%d (was %d)\n", nrow(disc_sens), nrow(disc)))
  cat("Groups:\n"); print(table(disc_sens$Group))
}

# ─── Helper: RF-LOOCV ────────────────────────────────────────────────────────
rf_loocv_auc <- function(df, ntree = 1000) {
  feat_cols <- setdiff(colnames(df), c("SampleID", "Group"))
  df_imp    <- df
  for (col in feat_cols) {
    df_imp[[col]][is.na(df_imp[[col]])] <- min(df_imp[[col]], na.rm=TRUE) / 2
  }
  X_raw <- df_imp[, feat_cols]
  y_fac <- factor(df_imp$Group, levels=c("Control","Infected"))
  cn    <- make.unique(gsub("[^A-Za-z0-9]","_", colnames(X_raw)), sep="_")
  cn    <- ifelse(grepl("^[0-9]",cn), paste0("F",cn), cn)
  colnames(X_raw) <- cn
  X_mat <- scale(log1p(as.matrix(X_raw)))
  n     <- nrow(X_mat)

  prob_vec <- numeric(n)
  for (i in seq_len(n)) {
    rf   <- randomForest(x = X_mat[-i,,drop=FALSE], y = y_fac[-i],
                          ntree = ntree, mtry = max(2, floor(sqrt(ncol(X_mat)))))
    prob_vec[i] <- predict(rf, newdata=X_mat[i,,drop=FALSE], type="prob")[,"Infected"]
  }
  roc_obj   <- roc(y_fac, prob_vec, levels=c("Control","Infected"),
                   direction="<", quiet=TRUE)
  ci_delong <- ci.auc(roc_obj, method="delong")
  list(auc=as.numeric(auc(roc_obj)), lo=ci_delong[1], hi=ci_delong[3],
       roc=roc_obj, n=n)
}

# ─── Primary analysis (all controls) ─────────────────────────────────────────
cat("Running primary RF-LOOCV (all controls, n=", nrow(disc), ")...\n")
res_primary <- rf_loocv_auc(disc)
cat(sprintf("  Primary AUC = %.3f (%.3f–%.3f)\n",
            res_primary$auc, res_primary$lo, res_primary$hi))

# ─── Sensitivity analysis (CAA+ excluded) ────────────────────────────────────
if (nrow(disc_sens) < nrow(disc)) {
  cat(sprintf("Running sensitivity RF-LOOCV (CAA+ excluded, n=%d)...\n", nrow(disc_sens)))
  res_sens <- rf_loocv_auc(disc_sens)
  cat(sprintf("  Sensitivity AUC = %.3f (%.3f–%.3f)\n",
              res_sens$auc, res_sens$lo, res_sens$hi))
} else {
  res_sens <- res_primary
  cat("Sensitivity = Primary (no samples excluded).\n")
}

# ─── Results table ───────────────────────────────────────────────────────────
sens_tbl <- tibble(
  Analysis = c("Primary (all controls)", "Sensitivity (CAA+ excluded)"),
  N        = c(res_primary$n, res_sens$n),
  AUC      = round(c(res_primary$auc, res_sens$auc), 3),
  CI_lo    = round(c(res_primary$lo,  res_sens$lo),  3),
  CI_hi    = round(c(res_primary$hi,  res_sens$hi),  3)
)
write_csv(sens_tbl, "tables/TableS_Sensitivity_CAA.csv")
cat("\n✓ TableS_Sensitivity_CAA.csv saved\n")
print(sens_tbl)

# ─── Overlay ROC plot ────────────────────────────────────────────────────────
roc_compare <- bind_rows(
  tibble(FPR=1-res_primary$roc$specificities,
         TPR=res_primary$roc$sensitivities,
         Analysis=sprintf("Primary (n=%d, AUC=%.3f)", res_primary$n, res_primary$auc)),
  tibble(FPR=1-res_sens$roc$specificities,
         TPR=res_sens$roc$sensitivities,
         Analysis=sprintf("CAA+ excluded (n=%d, AUC=%.3f)", res_sens$n, res_sens$auc))
)

p_sens <- ggplot(roc_compare, aes(x=FPR, y=TPR, color=Analysis, linetype=Analysis)) +
  geom_line(linewidth=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray60") +
  scale_color_manual(values=c("#2E86AB","#D85A30")) +
  labs(title    = "Sensitivity analysis — Discovery batch",
       subtitle = "7 CAA-positive KK− controls excluded vs. retained",
       x = "1 − Specificity", y = "Sensitivity", color=NULL, linetype=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",
        plot.title      = element_text(face="bold"),
        plot.subtitle   = element_text(size=9, color="gray40"))

ggsave("figures/FigS_Sensitivity_CAA.pdf", p_sens, width=6.5, height=5.5, dpi=300)
ggsave("figures/FigS_Sensitivity_CAA.png", p_sens, width=6.5, height=5.5, dpi=300)
cat("✓ FigS_Sensitivity_CAA saved\n")
