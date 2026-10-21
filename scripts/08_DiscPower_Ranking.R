##############################################################################
#  08_DiscPower_Ranking.R
#
#  Discriminatory power of each chemical class of differentially expressed
#  metabolites — bubble plot + pathway ranking bar chart.
#
#  For each chemical class:
#    - Individual AUC per metabolite (pROC, DeLong CI)
#    - Class-level mean AUC ± SE
#    - Bubble size = number of DEMs in class
#
#  Input  : data/DEM_annotation_reference.csv
#            data/Discovery_ML_input.csv
#            data/Validation_ML_input.csv
#            data/Polyparasitism_ML_input.csv
#  Output : figures/Fig_DiscPower_Bubble.pdf
#            figures/Fig_PathwayRanking_Bar.pdf
#            tables/Table_DiscPower_byClass.csv
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(ggrepel)
})

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ─── Load annotation + data ───────────────────────────────────────────────────
ann  <- read_csv("data/DEM_annotation_reference.csv", show_col_types = FALSE)
disc <- read_csv("data/Discovery_ML_input.csv",       show_col_types = FALSE)
val  <- read_csv("data/Validation_ML_input.csv",      show_col_types = FALSE)
poly <- read_csv("data/Polyparasitism_ML_input.csv",  show_col_types = FALSE)

batch_list <- list(Discovery = disc, Validation = val, Polyparasitism = poly)
for (nm in names(batch_list))
  batch_list[[nm]]$Group <- factor(batch_list[[nm]]$Group, levels = c("Control","Infected"))

# ─── Individual metabolite AUC per batch ─────────────────────────────────────
get_auc <- function(df, feat) {
  if (!feat %in% colnames(df)) return(NA_real_)
  vals <- df[[feat]]
  if (sum(!is.na(vals)) < 5) return(NA_real_)
  tryCatch({
    r <- roc(df$Group, vals, levels=c("Control","Infected"),
             direction="<", quiet=TRUE)
    as.numeric(auc(r))
  }, error=function(e) NA_real_)
}

feat_auc <- ann %>%
  rowwise() %>%
  mutate(
    AUC = get_auc(batch_list[[Batch]], Feature)
  ) %>%
  ungroup() %>%
  filter(!is.na(AUC))

cat("Feature-level AUC computed for n =", nrow(feat_auc), "metabolites\n")

# ─── Class-level summary ─────────────────────────────────────────────────────
class_summary <- feat_auc %>%
  group_by(ChemClass) %>%
  summarise(
    n_DEM   = n(),
    mean_AUC = mean(AUC, na.rm=TRUE),
    sd_AUC   = sd(AUC,   na.rm=TRUE),
    se_AUC   = sd_AUC / sqrt(n_DEM),
    n_TyrPath = if ("Is_Tyrosine_Pathway" %in% colnames(.))
                  sum(Is_Tyrosine_Pathway, na.rm=TRUE) else 0L,
    .groups = "drop"
  ) %>%
  mutate(
    pct_Tyr  = n_TyrPath / n_DEM * 100,
    Is_Tyr   = ChemClass %in% unique(ann$ChemClass[ann$Is_Tyrosine_Pathway == TRUE])
  ) %>%
  arrange(desc(mean_AUC))

write_csv(class_summary, "tables/Table_DiscPower_byClass.csv")
cat("✓ Table_DiscPower_byClass.csv saved\n")
print(class_summary)

# ─── Bubble plot: mean AUC vs #DEMs ──────────────────────────────────────────
p_bubble <- ggplot(class_summary,
                   aes(x = n_DEM, y = mean_AUC,
                       size = n_DEM, color = Is_Tyr,
                       label = ChemClass)) +
  geom_point(alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_AUC - se_AUC, ymax = mean_AUC + se_AUC),
                width = 0, linewidth = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  geom_text_repel(size = 3.2, max.overlaps = 15, box.padding = 0.4) +
  scale_size_continuous(range = c(3, 12), name = "# DEMs") +
  scale_color_manual(values = c("TRUE" = "#7F77DD", "FALSE" = "#B4B2A9"),
                     labels = c("TRUE" = "Tyrosine pathway", "FALSE" = "Other"),
                     name   = NULL) +
  scale_y_continuous(limits = c(0.45, 1.0),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title    = "Discriminatory power by chemical class",
    subtitle = "Mean AUC ± SE across individual DEMs per class | DeLong method",
    x = "Number of DEMs in class",
    y = "Mean AUC"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title   = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        legend.position = "bottom")

ggsave("figures/Fig_DiscPower_Bubble.pdf", p_bubble, width = 9, height = 6.5, dpi = 300)
ggsave("figures/Fig_DiscPower_Bubble.png", p_bubble, width = 9, height = 6.5, dpi = 300)
cat("✓ Fig_DiscPower_Bubble saved\n")

# ─── Pathway ranking bar chart ────────────────────────────────────────────────
p_bar <- class_summary %>%
  mutate(ChemClass = fct_reorder(ChemClass, mean_AUC)) %>%
  ggplot(aes(x = mean_AUC, y = ChemClass, fill = Is_Tyr)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(xmin = mean_AUC - se_AUC, xmax = mean_AUC + se_AUC),
                width = 0.3, linewidth = 0.7) +
  geom_text(aes(label = sprintf("n=%d  %.3f", n_DEM, mean_AUC)),
            hjust = -0.1, size = 3.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0.4, 1.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_fill_manual(values = c("TRUE" = "#7F77DD", "FALSE" = "#B4B2A9"),
                    labels = c("TRUE" = "Tyrosine pathway", "FALSE" = "Other"),
                    name   = NULL) +
  labs(
    title    = "Chemical class ranking by discriminatory power",
    subtitle = "Mean AUC ± SE | n = number of DEMs",
    x = "Mean AUC", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 9, color = "gray40"),
        legend.position = "bottom")

ggsave("figures/Fig_PathwayRanking_Bar.pdf", p_bar,
       width = 9, height = max(5, nrow(class_summary) * 0.6), dpi = 300)
ggsave("figures/Fig_PathwayRanking_Bar.png", p_bar,
       width = 9, height = max(5, nrow(class_summary) * 0.6), dpi = 300)
cat("✓ Fig_PathwayRanking_Bar saved\n")
