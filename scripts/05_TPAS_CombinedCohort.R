##############################################################################
#  05_TPAS_CombinedCohort.R
#
#  Combined-cohort LOOCV on pooled TPAS scores.
#  Uses batch-specific TPAS from script 03 — does NOT re-pool raw features
#  (which would cause direction inversion across batches).
#
#  Model: logistic regression (TPAS ~ Group) with LOOCV
#
#  Input  : data/TPAS_scores_allbatches.csv  (from script 03)
#  Output : figures/Fig4b_Combined_ROC.pdf
#
#  Expected: AUC = 0.847 (95% CI: 0.786–0.908)
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
})

dir.create("figures", showWarnings = FALSE)

# ─── Load TPAS scores ────────────────────────────────────────────────────────
tpas <- read_csv("data/TPAS_scores_allbatches.csv", show_col_types = FALSE)
tpas$Group <- factor(tpas$Group, levels = c("Control", "Infected"))
cat(sprintf("Loaded TPAS scores: n=%d | Batches: %s\n",
            nrow(tpas), paste(unique(tpas$Batch), collapse=", ")))
cat("Group distribution:\n"); print(table(tpas$Group, tpas$Batch))

# ─── Combined LOOCV (logistic regression) ────────────────────────────────────
n         <- nrow(tpas)
prob_loocv <- numeric(n)

for (i in seq_len(n)) {
  train    <- tpas[-i, ]
  test_row <- tpas[i, ]

  fit <- glm(Group ~ TPAS, data = train, family = binomial())
  prob_loocv[i] <- predict(fit, newdata = test_row, type = "response")
}

roc_obj  <- roc(tpas$Group, prob_loocv, levels = c("Control","Infected"),
                direction = "<", quiet = TRUE)
ci_delong <- ci.auc(roc_obj, method = "delong")
auc_val   <- as.numeric(auc(roc_obj))

# Optimal threshold (Youden index)
coords_best <- coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"),
                       best.method = "youden")
thresh <- coords_best["threshold"]
pred_class <- ifelse(prob_loocv >= as.numeric(thresh), "Infected", "Control")

sens <- mean(pred_class[tpas$Group == "Infected"] == "Infected")
spec <- mean(pred_class[tpas$Group == "Control"]  == "Control")

cat(sprintf("\n── Combined-cohort LOOCV results ──────────────────────────────\n"))
cat(sprintf("  AUC         = %.3f  (95%% CI: %.3f–%.3f)\n",
            auc_val, ci_delong[1], ci_delong[3]))
cat(sprintf("  Sensitivity = %.3f\n", sens))
cat(sprintf("  Specificity = %.3f\n", spec))
cat(sprintf("  Threshold   = %.3f (Youden index)\n\n", as.numeric(thresh)))

# ─── ROC plot ────────────────────────────────────────────────────────────────
roc_df <- tibble(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#7F77DD", linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(data = tibble(FPR = 1 - as.numeric(coords_best["specificity"]),
                            TPR = as.numeric(coords_best["sensitivity"])),
             aes(x = FPR, y = TPR), color = "#D85A30", size = 4, shape = 19) +
  annotate("text", x = 0.65, y = 0.25,
           label = sprintf("AUC = %.3f\n95%% CI: %.3f–%.3f\nSens=%.3f, Spec=%.3f",
                           auc_val, ci_delong[1], ci_delong[3], sens, spec),
           size = 3.8, hjust = 0, color = "#333333") +
  labs(
    title    = "Combined-cohort LOOCV — TPAS logistic model",
    subtitle = sprintf("n=%d pooled across 3 batches | logistic regression | Youden threshold = %.3f",
                       n, as.numeric(thresh)),
    x = "1 − Specificity (FPR)",
    y = "Sensitivity (TPR)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"))

ggsave("figures/Fig4b_Combined_ROC.pdf", p_roc, width = 6, height = 5.5, dpi = 300)
ggsave("figures/Fig4b_Combined_ROC.png", p_roc, width = 6, height = 5.5, dpi = 300)
cat("✓ Fig4b_Combined_ROC saved\n")
