##############################################################################
#  01_Boxplots_Tyrosine.R
#
#  Tyrosine pathway differentially expressed metabolites — boxplots
#  Discovery, Validation, and Polyparasitism batches
#
#  Input  : data/Discovery_ML_input.csv
#            data/Validation_ML_input.csv
#            data/Polyparasitism_ML_input.csv
#  Output : figures/Fig3_Boxplots_TyrPathway.pdf
#            figures/Fig3_Boxplots_TyrPathway.png
#
#  DEM selection thresholds:
#    Discovery       FC ≥ 1.2, p < 0.1  
#    Validation      FC ≥ 1.5, p < 0.05
#    Polyparasitism  FC ≥ 2.0, p < 0.05
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

col_infected <- "#D85A30"
col_control  <- "#1D9E75"

# ─── Helper: boxplot for one batch ──────────────────────────────────────────
make_boxplot <- function(df, tyr_cols, label_map, title_text, subtitle_text,
                          strip_fill = "#EEEDFE") {
  present <- tyr_cols[tyr_cols %in% colnames(df)]
  if (length(present) == 0) {
    warning("No tyrosine columns found for: ", title_text); return(NULL)
  }
  long <- df %>%
    select(SampleID, Group, all_of(present)) %>%
    pivot_longer(cols = all_of(present), names_to = "Feature", values_to = "Intensity") %>%
    mutate(
      Log_I = log1p(Intensity),
      Label = ifelse(Feature %in% names(label_map), label_map[Feature],
                     gsub("_", " ", Feature)),
      Label = factor(Label, levels = unname(label_map[present[present %in% names(label_map)]]))
    )

  stats <- long %>%
    group_by(Label) %>%
    summarise(
      p_val = tryCatch(
        wilcox.test(Log_I[Group == "Infected"], Log_I[Group == "Control"],
                    exact = FALSE)$p.value,
        error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(sig = case_when(
      p_val < 0.001 ~ "***", p_val < 0.01 ~ "**",
      p_val < 0.05  ~ "*",  TRUE ~ "ns")) %>%
    left_join(
      long %>% group_by(Label) %>%
        summarise(y_pos = max(Log_I, na.rm = TRUE) * 1.10, .groups = "drop"),
      by = "Label")

  ggplot(long, aes(x = Group, y = Log_I, fill = Group)) +
    geom_boxplot(width = 0.50, outlier.shape = NA, alpha = 0.72, linewidth = 0.5) +
    geom_jitter(aes(color = Group), width = 0.14, size = 2.2, alpha = 0.88) +
    geom_text(data = stats, aes(x = 1.5, y = y_pos, label = sig),
              size = 5, fontface = "bold", inherit.aes = FALSE) +
    facet_wrap(~Label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values  = c("Control" = col_control, "Infected" = col_infected), guide = "none") +
    scale_color_manual(values = c("Control" = col_control, "Infected" = col_infected), guide = "none") +
    labs(title = title_text, subtitle = subtitle_text,
         x = NULL, y = "log(intensity + 1)") +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 9,  color = "gray40"),
      strip.text       = element_text(size = 9,  face = "bold.italic"),
      strip.background = element_rect(fill = strip_fill, color = NA),
      axis.text.x      = element_text(size = 10, face = "bold")
    )
}

# ─── Load data ───────────────────────────────────────────────────────────────
disc <- read_csv("data/Discovery_ML_input.csv",     show_col_types = FALSE)
val  <- read_csv("data/Validation_ML_input.csv",    show_col_types = FALSE)
poly <- read_csv("data/Polyparasitism_ML_input.csv", show_col_types = FALSE)

for (d in list(disc, val, poly)) d$Group <- factor(d$Group, levels = c("Control", "Infected"))
disc$Group <- factor(disc$Group, levels = c("Control", "Infected"))
val$Group  <- factor(val$Group,  levels = c("Control", "Infected"))
poly$Group <- factor(poly$Group, levels = c("Control", "Infected"))

# ─── Discovery: FC≥1.2, p<0.1 ───────────────────────────────────────────────
disc_tyr <- c("L_Tyrosine", "L_DOPA", "Dopamine")
disc_labels <- c(
  "L_Tyrosine" = "L-Tyrosine\nFC=1.22↑ p=0.062",
  "L_DOPA"     = "L-DOPA\nFC=0.33↓ p=0.004",
  "Dopamine"   = "Dopamine\nFC=0.63↓ p=0.064"
)
fig_disc <- make_boxplot(disc, disc_tyr, disc_labels,
  "Tyrosine pathway — Discovery batch",
  "n=20 (10 infected, 10 controls) | FC≥1.2, p<0.1 | Wilcoxon rank-sum",
  strip_fill = "#E1F5EE")

# ─── Validation: FC≥1.5, p<0.05 ─────────────────────────────────────────────
val_tyr <- c("L_Tyrosine", "3_O_Methyldopamine", "trans_Cinnamic_acid", "Fumarate")
val_labels <- c(
  "L_Tyrosine"          = "L-Tyrosine\nFC=2.44↑ p<0.001",
  "3_O_Methyldopamine"  = "3-O-Methyldopamine\nFC=2.22↑ p<0.001",
  "trans_Cinnamic_acid" = "trans-Cinnamic acid\nFC=2.88↑ p=0.002",
  "Fumarate"            = "Fumarate\nFC=1.57↑ p=0.019"
)
fig_val <- make_boxplot(val, val_tyr, val_labels,
  "Tyrosine pathway — Validation batch",
  "n=37 (23 infected, 14 controls) | FC≥1.5, p<0.05 | Wilcoxon rank-sum",
  strip_fill = "#EEEDFE")

# ─── Polyparasitism: FC≥2.0, p<0.05 ─────────────────────────────────────────
poly_tyr <- c("Tyramine", "Hydroxyphenylacetylglycine",
              "X_4_Hydroxy_3_methoxyphenyl_ethanol",
              "4_Hydroxyphenylethanol", "4_Methylcatechol_1", "Gentisate_aldehyde")
poly_labels <- c(
  "Tyramine"                             = "Tyramine\nFC=2.24↑ p=0.006",
  "Hydroxyphenylacetylglycine"           = "Hydroxyphenyl-\nacetylglycine\nFC=3.08↑ p=0.016",
  "X_4_Hydroxy_3_methoxyphenyl_ethanol"  = "(4-OH-3-methoxy-\nphenyl)ethanol\nFC=2.84↑ p=0.003",
  "4_Hydroxyphenylethanol"               = "4-Hydroxyphenyl-\nethanol\nFC=2.46↑ p=0.012",
  "4_Methylcatechol_1"                   = "4-Methylcatechol\nFC=2.06↑ p=0.004",
  "Gentisate_aldehyde"                   = "Gentisate aldehyde\nFC=2.31↑ p=0.029"
)
fig_poly <- make_boxplot(poly, poly_tyr, poly_labels,
  "Tyrosine pathway — Polyparasitism batch",
  "n=12 (7 infected, 5 controls) | FC≥2.0, p<0.05 | Wilcoxon rank-sum",
  strip_fill = "#FAECE7")

# ─── Combine and save ────────────────────────────────────────────────────────
figs <- Filter(Negate(is.null), list(fig_disc, fig_val, fig_poly))
if (length(figs) > 0) {
  combined <- wrap_plots(figs, ncol = 1)
  ggsave("figures/Fig3_Boxplots_TyrPathway.pdf", combined,
         width = 14, height = 5.5 * length(figs), dpi = 300)
  ggsave("figures/Fig3_Boxplots_TyrPathway.png", combined,
         width = 14, height = 5.5 * length(figs), dpi = 300)
  cat("✓ Fig3_Boxplots_TyrPathway saved (", length(figs), "panels)\n")
}
