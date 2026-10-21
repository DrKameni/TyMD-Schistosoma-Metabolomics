##############################################################################
#  06_Clinical_Correlations.R
#
#  Spearman correlation of TPAS with clinical variables
#  (EPG, spleen score, liver score, age, sex) per batch and combined.
#  Multiple testing correction: Benjamini-Hochberg FDR.
#
#  Input  : data/TPAS_scores_allbatches.csv  (from script 03)
#            data/Clinical_master.csv
#  Output : figures/Fig_TPAS_Clinical_Heatmap.pdf
#            tables/Table_TPAS_Clinical_Correlations.csv
#
#  Expected: TPAS vs EPG rho ~ 0.41, p < 0.001 (all batches pooled)
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ─── Load data ────────────────────────────────────────────────────────────────
tpas <- read_csv("data/TPAS_scores_allbatches.csv", show_col_types = FALSE)
clin <- read_csv("data/Clinical_master.csv",         show_col_types = FALSE)

cat("TPAS n =", nrow(tpas), "| Clinical n =", nrow(clin), "\n")
cat("Clinical columns:", paste(colnames(clin), collapse=", "), "\n\n")

# ─── Join on SampleID + Batch ────────────────────────────────────────────────
# Must join on BOTH columns to avoid many-to-many when samples appear
# across batches in the clinical file.
merged <- left_join(tpas, clin, by = c("SampleID", "Batch"))
cat("After join: n =", nrow(merged), "\n")

# ─── Clinical variables to correlate ─────────────────────────────────────────
clin_vars <- intersect(
  c("EPG", "Spleen_score", "Liver_score", "Age", "Sex_numeric",
    "Worm_burden", "KK_count"),
  colnames(merged)
)
cat("Clinical variables found:", paste(clin_vars, collapse=", "), "\n\n")

if (length(clin_vars) == 0) stop("No clinical variables found after join. Check Clinical_master.csv column names.")

# ─── Spearman correlations (per batch + combined) ─────────────────────────────
batches <- c(unique(as.character(tpas$Batch)), "Combined")

corr_results <- map_dfr(batches, function(b) {
  df_sub <- if (b == "Combined") merged else filter(merged, Batch == b)

  map_dfr(clin_vars, function(cv) {
    vals <- df_sub[[cv]]
    tpas_v <- df_sub$TPAS
    ok    <- !is.na(vals) & !is.na(tpas_v)
    if (sum(ok) < 5) return(NULL)

    ct <- tryCatch(
      cor.test(tpas_v[ok], vals[ok], method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(ct)) return(NULL)

    tibble(
      Batch    = b,
      Variable = cv,
      n        = sum(ok),
      rho      = round(ct$estimate, 3),
      p_raw    = ct$p.value
    )
  })
})

# BH-FDR correction across all tests
corr_results <- corr_results %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH"),
         sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                           p_adj < 0.05  ~ "*",  TRUE ~ "ns"))

write_csv(corr_results, "tables/Table_TPAS_Clinical_Correlations.csv")
cat("✓ Table_TPAS_Clinical_Correlations.csv saved\n")
print(corr_results %>% filter(sig != "ns") %>% arrange(p_adj))

# ─── Correlation heatmap ─────────────────────────────────────────────────────
heat_df <- corr_results %>%
  mutate(
    Batch    = factor(Batch, levels = batches),
    Variable = factor(Variable, levels = clin_vars)
  )

p_heat <- ggplot(heat_df, aes(x = Batch, y = Variable, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), "\n", sig)),
            size = 3.2, lineheight = 0.9) +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#D85A30",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman ρ") +
  labs(
    title    = "TPAS correlation with clinical variables",
    subtitle = "Spearman ρ | BH-FDR corrected | ***p<0.001  **p<0.01  *p<0.05",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "gray40"),
    axis.text.x      = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.position  = "right"
  )

ggsave("figures/Fig_TPAS_Clinical_Heatmap.pdf", p_heat,
       width = 8, height = max(4, length(clin_vars) * 0.8), dpi = 300)
ggsave("figures/Fig_TPAS_Clinical_Heatmap.png", p_heat,
       width = 8, height = max(4, length(clin_vars) * 0.8), dpi = 300)
cat("✓ Fig_TPAS_Clinical_Heatmap saved\n")

# ─── Scatter: TPAS vs EPG (all batches) ──────────────────────────────────────
if ("EPG" %in% colnames(merged)) {
  merged_epg <- filter(merged, !is.na(EPG), EPG > 0)
  rho_epg  <- cor.test(merged_epg$TPAS, log(merged_epg$EPG),
                        method = "spearman", exact = FALSE)

  p_scatter <- ggplot(merged_epg, aes(x = log(EPG), y = TPAS, color = Batch)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "gray30", linetype = "dashed") +
    scale_color_manual(values = c("Discovery"="#2E86AB","Validation"="#1D9E75",
                                   "Polyparasitism"="#D85A30")) +
    annotate("text", x = min(log(merged_epg$EPG), na.rm=TRUE),
             y = max(merged_epg$TPAS, na.rm=TRUE),
             label = sprintf("ρ = %.2f\np = %.3g", rho_epg$estimate, rho_epg$p.value),
             hjust = 0, vjust = 1, size = 4, fontface = "bold") +
    labs(title    = "TPAS vs egg burden (EPG)",
         subtitle = "Spearman correlation | log-transformed EPG | KK+ samples only",
         x = "log(EPG)", y = "TPAS") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave("figures/Fig_TPAS_vs_EPG.pdf", p_scatter, width = 7, height = 5.5, dpi = 300)
  ggsave("figures/Fig_TPAS_vs_EPG.png", p_scatter, width = 7, height = 5.5, dpi = 300)
  cat("✓ Fig_TPAS_vs_EPG saved\n")
}
