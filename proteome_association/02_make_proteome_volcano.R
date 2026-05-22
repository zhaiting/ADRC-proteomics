#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    val <- sub("^--[^=]+=", "", arg)
    out[[gsub("-", "_", key)]] <- val
  }
  out
}

arg <- parse_args()

DREAM_OUT_DIR <- if (exists("DREAM_OUT_DIR", inherits = TRUE)) {
  get("DREAM_OUT_DIR", inherits = TRUE)
} else if (!is.null(arg$dream_dir)) {
  arg$dream_dir
} else {
  stop("Provide --dream-dir=<path> or set DREAM_OUT_DIR before sourcing.")
}

VOLCANO_PLOT_DIR <- if (exists("VOLCANO_PLOT_DIR", inherits = TRUE)) {
  get("VOLCANO_PLOT_DIR", inherits = TRUE)
} else if (!is.null(arg$plot_dir)) {
  arg$plot_dir
} else {
  file.path(DREAM_OUT_DIR, "volcano_plots")
}
dir.create(VOLCANO_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

PLOT_P_THRESHOLD <- if (exists("PLOT_P_THRESHOLD", inherits = TRUE)) {
  get("PLOT_P_THRESHOLD", inherits = TRUE)
} else if (!is.null(arg$p_threshold)) {
  as.numeric(arg$p_threshold)
} else {
  0.05
}

PLOT_LOG2FC_GUIDE <- if (exists("PLOT_LOG2FC_GUIDE", inherits = TRUE)) {
  get("PLOT_LOG2FC_GUIDE", inherits = TRUE)
} else if (!is.null(arg$log2fc_guide)) {
  as.numeric(arg$log2fc_guide)
} else {
  log2(1.1)
}

plot_color <- function(name, fallback) {
  if (exists(name, inherits = TRUE)) get(name, inherits = TRUE) else fallback
}

contrast_files <- list.files(
  DREAM_OUT_DIR,
  pattern = "^toptable_.*_vs_HC[.]csv$",
  full.names = FALSE
)
if (length(contrast_files) == 0L) {
  stop("No toptable_*_vs_HC.csv files found in: ", DREAM_OUT_DIR)
}

contrasts <- tibble(
  file = contrast_files,
  contrast = sub("^toptable_(.*)[.]csv$", "\\1", contrast_files)
)

format_title <- function(comparison_group, reference_group) {
  pretty_comparison <- gsub("_spectrum", "-spectrum", comparison_group)
  pretty_comparison <- gsub("_", "-", pretty_comparison)
  paste(pretty_comparison, "vs", reference_group)
}

read_contrast <- function(file, contrast) {
  df <- read_csv(file.path(DREAM_OUT_DIR, file), show_col_types = FALSE)
  comparison_group <- unique(df$comparison_group)[1]
  reference_group <- unique(df$reference_group)[1]
  n_col <- paste0("n_", comparison_group)

  df %>%
    mutate(
      contrast = contrast,
      title = format_title(comparison_group, reference_group),
      n_comparison = .data[[n_col]],
      neg_log10_p = -log10(pmax(P.Value, .Machine$double.xmin)),
      sig = adj.P.Val < 0.05,
      effect_bin = case_when(
        P.Value < PLOT_P_THRESHOLD & logFC >= PLOT_LOG2FC_GUIDE  ~ "Higher than cutoff",
        P.Value < PLOT_P_THRESHOLD & logFC <= -PLOT_LOG2FC_GUIDE ~ "Lower than cutoff",
        TRUE                                                     ~ "Not highlighted"
      ),
      effect_bin = factor(
        effect_bin,
        levels = c("Higher than cutoff", "Lower than cutoff", "Not highlighted")
      ),
      label = case_when(
        !is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "" ~ EntrezGeneSymbol,
        !is.na(Target) & Target != "" ~ Target,
        TRUE ~ SomaKey
      )
    )
}

all_res <- purrr::pmap_dfr(contrasts, read_contrast)

summary_tbl <- all_res %>%
  group_by(contrast, title, comparison_group, reference_group) %>%
  summarise(
    n_proteins = n(),
    n_comparison = first(n_comparison),
    n_HC = first(n_HC),
    group_mode = first(group_mode),
    min_baseline_age = first(min_baseline_age),
    model_formula = first(model_formula),
    fdr_0_05 = sum(sig, na.rm = TRUE),
    higher = sum(sig & logFC > 0, na.rm = TRUE),
    lower = sum(sig & logFC < 0, na.rm = TRUE),
    min_adj_p = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(VOLCANO_PLOT_DIR, "volcano_summary.csv"))

sig_tbl <- all_res %>%
  filter(sig) %>%
  arrange(contrast, adj.P.Val) %>%
  select(contrast, comparison_group, reference_group, SomaKey, Target,
         TargetFullName, UniProt, EntrezGeneSymbol, logFC, P.Value, adj.P.Val,
         n_comparison, n_HC, model_formula)
write_csv(sig_tbl, file.path(VOLCANO_PLOT_DIR, "fdr_significant_proteins.csv"))

palette <- c(
  "Higher than cutoff" = plot_color("PLOT_COLOR_HIGHER", "#D55E00"),
  "Lower than cutoff" = plot_color("PLOT_COLOR_LOWER", "#0072B2"),
  "Not highlighted" = plot_color("PLOT_COLOR_NEUTRAL", "#a6adb7")
)

make_plot <- function(df) {
  label_df <- df %>%
    arrange(P.Value) %>%
    mutate(label_me = sig | row_number() <= 5) %>%
    filter(label_me)

  ggplot(df, aes(x = logFC, y = neg_log10_p)) +
    geom_hline(yintercept = -log10(PLOT_P_THRESHOLD), color = "grey75",
               linewidth = 0.35, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey78", linewidth = 0.35) +
    geom_vline(xintercept = c(-PLOT_LOG2FC_GUIDE, PLOT_LOG2FC_GUIDE),
               color = "grey65", linewidth = 0.3, linetype = "dotted") +
    geom_point(aes(color = effect_bin), alpha = 0.72, size = 1.45) +
    geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 3.2,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.35,
      seed = 1
    ) +
    scale_color_manual(values = palette, breaks = names(palette),
                       drop = FALSE, guide = "none") +
    labs(
      title = unique(df$title),
      x = "Model log2 fold change vs HC",
      y = "-log10 P value",
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

plots <- all_res %>%
  split(.$contrast) %>%
  lapply(make_plot)

for (contrast in names(plots)) {
  ggsave(
    file.path(VOLCANO_PLOT_DIR, sprintf("volcano_%s.png", contrast)),
    plots[[contrast]], width = 7.4, height = 5.6, dpi = 300
  )
  ggsave(
    file.path(VOLCANO_PLOT_DIR, sprintf("volcano_%s.pdf", contrast)),
    plots[[contrast]], width = 7.4, height = 5.6
  )
}

message("Wrote full-proteome volcano plots to: ", VOLCANO_PLOT_DIR)
