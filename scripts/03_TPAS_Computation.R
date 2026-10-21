##############################################################################
#  03_TPAS_Computation.R
#
#  Tyrosine Pathway Activity Score (TPAS) — direction-corrected composite
#  z-score computed independently within each batch.
#
#  Method:
#    1. For each DEM, compute within-batch z-score
#    2. Flip sign for metabolites with DOWN direction (FC < 1 in infected)
#       so that ALL features contribute positively when infection is present
#    3. TPAS = mean of direction-corrected z-scores across DEM features
#    4. Wilcoxon test (KK+ vs KK−) per batch + combined
#
#  Features per batch (shared across ≥2 batches define the core TyMD signature):
#    Discovery       : L_Tyrosine (↑), L_DOPA (↓), Dopamine (↓)
#    Validation      : L_Tyrosine (↑), 3_O_Methyldopamine (↑),
#                      trans_Cinnamic_acid (↑), Fumarate (↑)
#    Polyparasitism  : Tyramine (↑), Hydroxyphenylacetylglycine (↑),
#                      X_4_Hydroxy_3_methoxyphenyl_ethanol (↑),
#                      4_Hydroxyphenylethanol (↑), 4_Methylcatechol_1 (↑),
#                      Gentisate_aldehyde (↑)
#
#  Input  : data/Discovery_ML_input.csv
#            data/Validation_ML_input.csv
#            data/Polyparasitism_ML_input.csv
#  Output : data/TPAS_scores_allbatches.csv   ← used by scripts 04 and 05
#            figures/Fig6_TPAS_scores.pdf
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(coin)        # Wilcoxon exact permutation test (optional fallback)
})

dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

col_infected <- "#D85A30"
col_control  <- "#1D9E75"

# ─── Batch feature definitions ───────────────────────────────────────────────
# direction: +1 = up in infected, -1 = down in infected
batch_defs <- list(
  Discovery = list(
    file     = "data/Discovery_ML_input.csv",
    features = c("L_Tyrosine", "L_DOPA", "Dopamine"),
    dirs     = c(+1, -1, -1)          # L_DOPA and Dopamine DOWN in infected
  ),
  Validation = list(
    file     = "data/Validation_ML_input.csv",
    features = c("L_Tyrosine", "3_O_Methyldopamine",
                 "trans_Cinnamic_acid", "Fumarate"),
    dirs     = c(+1, +1, +1, +1)
  ),
  Polyparasitism = list(
    file     = "data/Polyparasitism_ML_input.csv",
    features = c("Tyramine", "Hydroxyphenylacetylglycine",
                 "X_4_Hydroxy_3_methoxyphenyl_ethanol",
                 "4_Hydroxyphenylethanol", "4_Methylcatechol_1",
                 "Gentisate_aldehyde"),
    dirs     = c(+1, +1, +1, +1, +1, +1)
  )
)

# ─── Compute TPAS per batch ──────────────────────────────────────────────────
tpas_all <- list()

for (batch_nm in names(batch_defs)) {
  def  <- batch_defs[[batch_nm]]
  df   <- read_csv(def$file, show_col_types = FALSE)
  df$Group <- factor(df$Group, levels = c("Control", "Infected"))

  # Keep only features present in this batch
  feats   <- def$features
  dirs    <- def$dirs
  present <- feats %in% colnames(df)
  if (!all(present)) {
    warning(sprintf("[%s] Missing features: %s",
                    batch_nm, paste(feats[!present], collapse = ", ")))
  }
  feats <- feats[present]
  dirs  <- dirs[present]

  # Missing value imputation: half-minimum per column
  for (col in feats) {
    df[[col]][is.na(df[[col]])] <- min(df[[col]], na.rm = TRUE) / 2
  }

  # Z-score within batch
  feat_mat <- as.matrix(df[, feats, drop = FALSE])
  feat_z   <- scale(log1p(feat_mat))          # z-score of log-transformed values

  # Direction correction
  for (j in seq_along(feats)) feat_z[, j] <- feat_z[, j] * dirs[j]

  # TPAS = row mean of direction-corrected z-scores
  tpas_vec <- rowMeans(feat_z, na.rm = TRUE)

  tpas_all[[batch_nm]] <- tibble(
    SampleID = df$SampleID,
    Group    = df$Group,
    Batch    = batch_nm,
    TPAS     = tpas_vec
  )

  # Wilcoxon test
  wt <- wilcox.test(TPAS ~ Group, data = tpas_all[[batch_nm]], exact = FALSE)
  cat(sprintf("[%s] TPAS Wilcoxon p = %.4f | Infected median = %.3f | Control median = %.3f\n",
              batch_nm,
              wt$p.value,
              median(tpas_vec[df$Group == "Infected"]),
              median(tpas_vec[df$Group == "Control"])))
}

# ─── Save TPAS scores (used by scripts 04 and 05) ───────────────────────────
tpas_combined <- bind_rows(tpas_all)
write_csv(tpas_combined, "data/TPAS_scores_allbatches.csv")
cat("\n✓ data/TPAS_scores_allbatches.csv saved\n")

# Combined Wilcoxon
wt_comb <- wilcox.test(TPAS ~ Group, data = tpas_combined, exact = FALSE)
cat(sprintf("Combined TPAS Wilcoxon p = %.4g\n\n", wt_comb$p.value))

# ─── TPAS boxplot per batch + combined ──────────────────────────────────────
plot_tpas <- function(df_sub, title_text, strip_fill = "#EEEDFE") {
  wt <- wilcox.test(TPAS ~ Group, data = df_sub, exact = FALSE)
  sig_label <- if (wt$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", wt$p.value)
  y_max <- max(df_sub$TPAS, na.rm = TRUE)

  ggplot(df_sub, aes(x = Group, y = TPAS, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.72, linewidth = 0.5) +
    geom_jitter(aes(color = Group), width = 0.15, size = 2.5, alpha = 0.85) +
    annotate("text", x = 1.5, y = y_max * 1.08,
             label = sig_label, size = 4.5, fontface = "bold") +
    scale_fill_manual(values  = c("Control" = col_control, "Infected" = col_infected), guide = "none") +
    scale_color_manual(values = c("Control" = col_control, "Infected" = col_infected), guide = "none") +
    labs(title = title_text, x = NULL, y = "TPAS (direction-corrected z-score)") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          axis.text.x   = element_text(size = 10, face = "bold"))
}

batch_fills <- c(Discovery = "#E1F5EE", Validation = "#EEEDFE", Polyparasitism = "#FAECE7")

figs_batch <- imap(tpas_all, function(df, nm) {
  plot_tpas(df,
    sprintf("%s (n=%d)", nm, nrow(df)),
    strip_fill = batch_fills[[nm]])
})

fig_combined <- plot_tpas(tpas_combined,
  sprintf("All batches combined (n=%d)", nrow(tpas_combined)))

fig_full <- wrap_plots(c(figs_batch, list(fig_combined)), nrow = 1)

ggsave("figures/Fig6_TPAS_scores.pdf", fig_full, width = 16, height = 5.5, dpi = 300)
ggsave("figures/Fig6_TPAS_scores.png", fig_full, width = 16, height = 5.5, dpi = 300)
cat("✓ Fig6_TPAS_scores saved\n")
