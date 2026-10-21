##############################################################################
#  04_MetaAnalysis_AUC.R
#
#  Random-effects meta-analysis of RF-LOOCV AUC across 3 batches
#  Method: DerSimonian-Laird (DL) with logit-AUC transformation
#
#  Input  : tables/Table2_RF_Performance.csv  (from script 02)
#  Output : figures/Fig7_MetaAnalysis_Forest.pdf
#            figures/Fig7_MetaAnalysis_Forest.png
#
#  Expected results:
#    Pooled AUC = 0.911 (95% CI: 0.823–0.960)
#    I²         = 0.0%
#    Q-test p   = 0.847
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(meta)
})

dir.create("figures", showWarnings = FALSE)

# ─── Load RF performance table ───────────────────────────────────────────────
perf <- read_csv("tables/Table2_RF_Performance.csv", show_col_types = FALSE)
cat("Loaded performance table:\n"); print(perf)

# ─── Logit transformation for meta-analysis ──────────────────────────────────
# AUC is bounded [0,1] → logit transform before DL pooling
logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Approximate SE of logit-AUC using CI endpoints (DeLong method)
perf_meta <- perf %>%
  mutate(
    logit_auc   = logit(AUC),
    logit_lo    = logit(pmax(CI_lo, 0.501)),   # avoid logit(0.5)=0 boundary
    logit_hi    = logit(pmin(CI_hi, 0.999)),
    se_logit    = (logit_hi - logit_lo) / (2 * 1.96)
  )

# DerSimonian-Laird meta-analysis on logit scale
ma <- metagen(
  TE    = perf_meta$logit_auc,
  seTE  = perf_meta$se_logit,
  studlab = perf_meta$Batch,
  sm    = "AUC",
  method.tau = "DL",
  random = TRUE,
  fixed  = FALSE,
  title  = "TyMD signature — RF-LOOCV AUC meta-analysis"
)

# Back-transform pooled estimate
pooled_auc <- inv_logit(ma$TE.random)
pooled_lo  <- inv_logit(ma$lower.random)
pooled_hi  <- inv_logit(ma$upper.random)

cat(sprintf("\n── Meta-analysis results ──────────────────────────────────────\n"))
cat(sprintf("  Pooled AUC = %.3f  (95%% CI: %.3f–%.3f)\n", pooled_auc, pooled_lo, pooled_hi))
cat(sprintf("  I²         = %.1f%%\n", ma$I2 * 100))
cat(sprintf("  Q-test p   = %.3f\n", ma$pval.Q))
cat(sprintf("  tau²       = %.4f\n", ma$tau2))
cat(sprintf("  k          = %d batches\n\n", ma$k))

# ─── ggplot2 forest plot ─────────────────────────────────────────────────────
batch_colors <- c("Discovery" = "#2E86AB", "Validation" = "#1D9E75",
                  "Polyparasitism" = "#D85A30", "Pooled" = "#333333")

# Study-level AUC back-transformed from logit
study_df <- perf_meta %>%
  mutate(
    AUC_bt   = inv_logit(logit_auc),
    CI_lo_bt = inv_logit(logit_auc - 1.96 * se_logit),
    CI_hi_bt = inv_logit(logit_auc + 1.96 * se_logit),
    Type     = "Study"
  ) %>%
  select(Batch, AUC_bt, CI_lo_bt, CI_hi_bt, N, Type)

# Pooled row
pooled_row <- tibble(
  Batch    = "Pooled (DL)",
  AUC_bt   = pooled_auc,
  CI_lo_bt = pooled_lo,
  CI_hi_bt = pooled_hi,
  N        = sum(perf$N),
  Type     = "Pooled"
)

forest_df <- bind_rows(study_df, pooled_row) %>%
  mutate(
    Batch = factor(Batch, levels = rev(c("Discovery","Validation",
                                          "Polyparasitism","Pooled (DL)"))),
    label = sprintf("%.3f (%.3f–%.3f)", AUC_bt, CI_lo_bt, CI_hi_bt)
  )

fp <- ggplot(forest_df,
             aes(x = AUC_bt, y = Batch,
                 xmin = CI_lo_bt, xmax = CI_hi_bt,
                 color = ifelse(Type == "Pooled", "Pooled", Batch),
                 shape = Type, size = Type)) +
  geom_vline(xintercept = 0.5,  linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = pooled_auc, linetype = "dotted", color = "#333333", linewidth = 0.8) +
  geom_errorbar(aes(xmin = CI_lo_bt, xmax = CI_hi_bt, y = Batch),
                height = 0.25, linewidth = 0.8) +
  geom_point() +
  geom_text(aes(x = 1.02, label = label), hjust = 0, size = 3.2, color = "black") +
  geom_text(aes(x = 0.45,  label = paste0("n=", N)),  hjust = 1, size = 3.2, color = "gray40") +
  scale_x_continuous(limits = c(0.4, 1.18), breaks = seq(0.4, 1.0, 0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = batch_colors, guide = "none") +
  scale_shape_manual(values = c("Study" = 15, "Pooled" = 18), guide = "none") +
  scale_size_manual(values  = c("Study" = 3,  "Pooled" = 5),  guide = "none") +
  annotate("text", x = 0.475, y = 0.55,
           label = sprintf("I²=%.1f%%  Q-p=%.3f", ma$I2 * 100, ma$pval.Q),
           size = 3.5, color = "gray30", hjust = 0) +
  labs(
    title    = "Meta-analysis of RF-LOOCV AUC — TyMD signature",
    subtitle = sprintf("Random effects (DerSimonian-Laird) | k=3 batches | Pooled AUC=%.3f (95%%CI:%.3f–%.3f)",
                       pooled_auc, pooled_lo, pooled_hi),
    x = "AUC (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y   = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)
  )

ggsave("figures/Fig7_MetaAnalysis_Forest.pdf", fp, width = 9, height = 5, dpi = 300)
ggsave("figures/Fig7_MetaAnalysis_Forest.png", fp, width = 9, height = 5, dpi = 300)
cat("✓ Fig7_MetaAnalysis_Forest saved\n")
