##############################################################################
#  07_ChemClass_Waffle.R
#
#  Chemical class distribution of differentially expressed metabolites (DEM)
#  — waffle charts per batch and combined.
#
#  Input  : data/DEM_annotation_reference.csv
#  Output : figures/FigS_ChemClass_Waffle.pdf
#
#  DEM_annotation_reference.csv must contain columns:
#    Feature, Batch, ChemClass, Is_Tyrosine_Pathway (TRUE/FALSE)
##############################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(waffle)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

# ─── Load annotation ─────────────────────────────────────────────────────────
ann <- read_csv("data/DEM_annotation_reference.csv", show_col_types = FALSE)

# Check required columns
required_cols <- c("Feature","Batch","ChemClass")
missing_cols  <- setdiff(required_cols, colnames(ann))
if (length(missing_cols) > 0) stop("Missing columns in DEM_annotation_reference.csv: ",
                                    paste(missing_cols, collapse=", "))

cat("DEM annotation loaded: n =", nrow(ann), "\n")
cat("Batches:", paste(unique(ann$Batch), collapse=", "), "\n")
cat("Chemical classes:", paste(unique(ann$ChemClass), collapse=", "), "\n\n")

# ─── Colour palette for chemical classes ─────────────────────────────────────
class_colors <- c(
  "Amino acids & peptides" = "#7F77DD",
  "Acylcarnitines"         = "#2E86AB",
  "Fatty acids"            = "#F4A261",
  "Organic acids"          = "#1D9E75",
  "Phenolic compounds"     = "#D85A30",
  "Nucleosides & purines"  = "#E9C46A",
  "Lipids"                 = "#264653",
  "Unknowns"               = "#B4B2A9",
  "Other"                  = "#CCD5AE"
)

# ─── Build waffle per batch ───────────────────────────────────────────────────
make_waffle <- function(df_sub, batch_nm) {
  counts <- df_sub %>%
    count(ChemClass, name = "n") %>%
    arrange(desc(n))

  vals <- setNames(counts$n, counts$ChemClass)
  cols <- class_colors[names(vals)]
  cols[is.na(cols)] <- "gray70"   # fallback for unlisted classes

  waffle(vals,
         rows         = 5,
         colors       = cols,
         title        = sprintf("%s (n=%d DEMs)", batch_nm, sum(vals)),
         xlab         = "1 square = 1 metabolite",
         legend_pos   = "bottom",
         size         = 0.4)
}

batch_plots <- map(c("Discovery","Validation","Polyparasitism"), function(b) {
  df_b <- filter(ann, Batch == b)
  if (nrow(df_b) == 0) return(NULL)
  make_waffle(df_b, b)
})

# Combined
combined_plot <- make_waffle(ann, "All batches combined")

# ─── Arrange and save ────────────────────────────────────────────────────────
all_plots <- Filter(Negate(is.null), batch_plots)
all_plots[[length(all_plots) + 1]] <- combined_plot

fig_waffle <- wrap_plots(all_plots, ncol = 2) +
  plot_annotation(
    title    = "Chemical class distribution of differentially expressed metabolites",
    subtitle = "Untargeted LC-MS/MS | ZIC-HILIC | Orbitrap",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  )

ggsave("figures/FigS_ChemClass_Waffle.pdf", fig_waffle,
       width = 14, height = 10, dpi = 300)
ggsave("figures/FigS_ChemClass_Waffle.png", fig_waffle,
       width = 14, height = 10, dpi = 300)
cat("✓ FigS_ChemClass_Waffle saved\n")

# ─── Summary table ───────────────────────────────────────────────────────────
summary_tbl <- ann %>%
  count(Batch, ChemClass, name = "n_DEM") %>%
  pivot_wider(names_from = Batch, values_from = n_DEM, values_fill = 0) %>%
  arrange(desc(rowSums(across(where(is.numeric)))))
print(summary_tbl)
