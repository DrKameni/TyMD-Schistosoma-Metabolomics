##############################################################################
#  02_ML_RandomForest_LOOCV.R
#
#  Random Forest Leave-One-Out Cross-Validation (RF-LOOCV) + SHAP importance
#  for each of the three batches independently.
#
#  Input  : data/Discovery_ML_input.csv
#            data/Validation_ML_input.csv
#            data/Polyparasitism_ML_input.csv
#  Output : figures/Fig4_RFLOOCV_ROC.pdf
#            figures/Fig5_SHAP_Importance.pdf
#            tables/Table2_RF_Performance.csv
#            tables/Table3_SHAP_Top15.csv
#
#  ⚠ LEAKAGE NOTE:
#  Volcano-based feature selection was performed on the FULL batch before LOOCV.
#  AUC values are therefore optimistic upper bounds, not unbiased estimates.
#  True generalisation is better captured by the meta-analysis (script 04)
#  or combined LOOCV (script 05). This is acknowledged in the manuscript.
#
#  Method:
#    - Missing value imputation: half-minimum per column (within batch)
#    - Log1p transformation + z-score scaling
#    - RF ntree=1000, mtry tuned by 3-fold CV
#    - SHAP values via iml::Predictor + iml::Shapley (sample-level)
#    - AUC 95% CI: DeLong method (pROC)
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(pROC)
  library(iml)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ─── Helper: half-min imputation (column-wise) ───────────────────────────────
impute_halfmin <- function(df) {
  num_cols <- names(df)[sapply(df, is.numeric)]
  for (col in num_cols) {
    df[[col]][is.na(df[[col]])] <- min(df[[col]], na.rm = TRUE) / 2
  }
  df
}

# ─── Helper: RF-LOOCV for one batch ─────────────────────────────────────────
run_rf_loocv <- function(df, batch_name, ntree = 1000) {
  cat(sprintf("\n=== %s RF-LOOCV ===\n", batch_name))

  feat_cols <- setdiff(colnames(df), c("SampleID", "Group"))
  df_imp    <- impute_halfmin(df)

  X_raw  <- df_imp[, feat_cols]
  y_fac  <- factor(df_imp$Group, levels = c("Control", "Infected"))

  # Clean column names for randomForest
  cn <- colnames(X_raw)
  cn <- gsub("[^A-Za-z0-9]", "_", cn)
  cn <- gsub("_{2,}", "_", cn)
  cn <- gsub("^_|_$", "", cn)
  cn <- ifelse(grepl("^[0-9]", cn), paste0("F", cn), cn)
  cn <- make.unique(cn, sep = "_")
  colnames(X_raw) <- cn

  X_mat <- scale(log1p(as.matrix(X_raw)))
  n     <- nrow(X_mat)
  cat(sprintf("  Samples: %d | Features: %d\n", n, ncol(X_mat)))

  # Tune mtry by 3-fold CV
  folds    <- sample(rep(1:3, length.out = n))
  mtry_try <- unique(c(2, 4, 6, 8, floor(sqrt(ncol(X_mat))), min(15, ncol(X_mat))))
  mtry_auc <- sapply(mtry_try, function(m) {
    mean(sapply(1:3, function(k) {
      rf  <- randomForest(x = X_mat[folds != k, , drop = FALSE],
                          y = y_fac[folds != k], ntree = 300, mtry = m)
      pr  <- predict(rf, newdata = X_mat[folds == k, , drop = FALSE], type = "prob")[, "Infected"]
      if (length(unique(y_fac[folds == k])) == 2)
        as.numeric(auc(roc(y_fac[folds == k], pr, levels = c("Control","Infected"),
                           direction = "<", quiet = TRUE)))
      else NA_real_
    }), na.rm = TRUE)
  })
  best_mtry <- mtry_try[which.max(mtry_auc)]
  cat(sprintf("  Best mtry: %d\n", best_mtry))

  # LOOCV
  prob_vec <- numeric(n)
  pred_vec <- character(n)
  for (i in seq_len(n)) {
    rf_i <- randomForest(x = X_mat[-i, , drop = FALSE],
                         y = y_fac[-i],
                         ntree = ntree, mtry = best_mtry)
    prob_vec[i] <- predict(rf_i, newdata = X_mat[i, , drop = FALSE], type = "prob")[, "Infected"]
    pred_vec[i] <- as.character(predict(rf_i, newdata = X_mat[i, , drop = FALSE]))
  }

  roc_obj <- roc(y_fac, prob_vec, levels = c("Control","Infected"),
                 direction = "<", quiet = TRUE)
  ci_obj  <- ci.auc(roc_obj, method = "delong")
  auc_val <- as.numeric(auc(roc_obj))
  sens    <- mean(pred_vec[y_fac == "Infected"]   == "Infected")
  spec    <- mean(pred_vec[y_fac == "Control"]    == "Control")

  cat(sprintf("  AUC = %.3f (%.3f–%.3f) | Sens = %.3f | Spec = %.3f\n",
              auc_val, ci_obj[1], ci_obj[3], sens, spec))

  # SHAP (on full model trained on all samples)
  rf_full   <- randomForest(x = X_mat, y = y_fac, ntree = ntree, mtry = best_mtry)
  pred_fn   <- function(model, newdata) predict(model, newdata = newdata, type = "prob")[, "Infected"]
  predictor <- Predictor$new(rf_full, data = as.data.frame(X_mat),
                             y = as.numeric(y_fac == "Infected"),
                             predict.function = pred_fn)

  shap_list <- lapply(seq_len(min(n, 30)), function(i) {
    tryCatch(Shapley$new(predictor, x.interest = as.data.frame(X_mat)[i, ])$results,
             error = function(e) NULL)
  })
  shap_df <- bind_rows(Filter(Negate(is.null), shap_list)) %>%
    group_by(feature) %>%
    summarise(mean_abs_shap = mean(abs(phi), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abs_shap)) %>%
    mutate(Batch = batch_name)

  list(
    batch    = batch_name,
    roc      = roc_obj,
    auc      = auc_val,
    ci_lo    = ci_obj[1],
    ci_hi    = ci_obj[3],
    sens     = sens,
    spec     = spec,
    shap     = shap_df,
    prob_vec = prob_vec,
    obs      = y_fac
  )
}

# ─── Load data ────────────────────────────────────────────────────────────────
disc <- read_csv("data/Discovery_ML_input.csv",      show_col_types = FALSE)
val  <- read_csv("data/Validation_ML_input.csv",     show_col_types = FALSE)
poly <- read_csv("data/Polyparasitism_ML_input.csv", show_col_types = FALSE)

# ─── Run LOOCV ────────────────────────────────────────────────────────────────
res_disc <- run_rf_loocv(disc, "Discovery")
res_val  <- run_rf_loocv(val,  "Validation")
res_poly <- run_rf_loocv(poly, "Polyparasitism")

# ─── Performance table ────────────────────────────────────────────────────────
perf_tbl <- map_dfr(list(res_disc, res_val, res_poly), function(r) {
  tibble(Batch       = r$batch,
         N           = length(r$obs),
         N_infected  = sum(r$obs == "Infected"),
         N_control   = sum(r$obs == "Control"),
         AUC         = round(r$auc, 3),
         CI_lo       = round(r$ci_lo, 3),
         CI_hi       = round(r$ci_hi, 3),
         Sensitivity = round(r$sens, 3),
         Specificity = round(r$spec, 3))
})
write_csv(perf_tbl, "tables/Table2_RF_Performance.csv")
cat("\n✓ Table2_RF_Performance.csv saved\n")
print(perf_tbl)

# ─── SHAP top-15 table ───────────────────────────────────────────────────────
shap_all <- bind_rows(res_disc$shap, res_val$shap, res_poly$shap) %>%
  group_by(Batch) %>%
  slice_max(mean_abs_shap, n = 15) %>%
  ungroup()
write_csv(shap_all, "tables/Table3_SHAP_Top15.csv")
cat("✓ Table3_SHAP_Top15.csv saved\n")

# ─── ROC plot ────────────────────────────────────────────────────────────────
batch_colors <- c("Discovery" = "#2E86AB", "Validation" = "#1D9E75", "Polyparasitism" = "#D85A30")

roc_plot <- local({
  roc_data <- imap_dfr(
    list(Discovery = res_disc, Validation = res_val, Polyparasitism = res_poly),
    function(r, nm) {
      tibble(Specificity = r$roc$specificities,
             Sensitivity = r$roc$sensitivities,
             Batch       = sprintf("%s (AUC=%.3f, 95%%CI:%.3f–%.3f)",
                                   nm, r$auc, r$ci_lo, r$ci_hi))
    })

  ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Batch)) +
    geom_line(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    scale_color_manual(values = setNames(
      batch_colors, sprintf("%s (AUC=%.3f, 95%%CI:%.3f–%.3f)",
                            names(batch_colors),
                            c(res_disc$auc, res_val$auc, res_poly$auc),
                            c(res_disc$ci_lo, res_val$ci_lo, res_poly$ci_lo),
                            c(res_disc$ci_hi, res_val$ci_hi, res_poly$ci_hi)))) +
    labs(title    = "RF-LOOCV ROC curves — TyMD signature",
         subtitle = "Three independent cohorts | DeLong 95% CI",
         x = "1 − Specificity", y = "Sensitivity", color = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 9),
          plot.title      = element_text(face = "bold"))
})

ggsave("figures/Fig4_RFLOOCV_ROC.pdf", roc_plot, width = 7, height = 6, dpi = 300)
ggsave("figures/Fig4_RFLOOCV_ROC.png", roc_plot, width = 7, height = 6, dpi = 300)
cat("✓ Fig4_RFLOOCV_ROC saved\n")

# ─── SHAP importance plot ────────────────────────────────────────────────────
shap_plot <- shap_all %>%
  group_by(Batch) %>%
  mutate(rank = row_number(),
         feature_label = gsub("_", " ", feature)) %>%
  ggplot(aes(x = mean_abs_shap, y = reorder(feature_label, mean_abs_shap),
             fill = Batch)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Batch, scales = "free_y") +
  scale_fill_manual(values = batch_colors) +
  labs(title    = "SHAP feature importance — top 15 per batch",
       subtitle = "Mean |SHAP| across samples; RF trained on all samples",
       x = "Mean |SHAP value|", y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold"),
        axis.text.y  = element_text(size = 8))

ggsave("figures/Fig5_SHAP_Importance.pdf", shap_plot, width = 14, height = 7, dpi = 300)
ggsave("figures/Fig5_SHAP_Importance.png", shap_plot, width = 14, height = 7, dpi = 300)
cat("✓ Fig5_SHAP_Importance saved\n")
