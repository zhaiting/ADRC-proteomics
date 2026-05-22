#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
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

PANEL_SOURCE <- if (exists("PANEL_SOURCE", inherits = TRUE)) {
  get("PANEL_SOURCE", inherits = TRUE)
} else if (!is.null(arg$panel_source)) {
  arg$panel_source
} else {
  stop("Provide --panel-source=<path> or set PANEL_SOURCE before sourcing.")
}

PANEL_OUT_DIR <- if (exists("PANEL_OUT_DIR", inherits = TRUE)) {
  get("PANEL_OUT_DIR", inherits = TRUE)
} else {
  file.path(DREAM_OUT_DIR, "panel_subsets")
}

PANEL_VOLCANO_PLOT_DIR <- if (exists("PANEL_VOLCANO_PLOT_DIR", inherits = TRUE)) {
  get("PANEL_VOLCANO_PLOT_DIR", inherits = TRUE)
} else {
  file.path(PANEL_OUT_DIR, "volcano_plots")
}
dir.create(PANEL_VOLCANO_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

PLOT_P_THRESHOLD <- if (exists("PLOT_P_THRESHOLD", inherits = TRUE)) {
  get("PLOT_P_THRESHOLD", inherits = TRUE)
} else {
  0.05
}

PLOT_LOG2FC_GUIDE <- if (exists("PLOT_LOG2FC_GUIDE", inherits = TRUE)) {
  get("PLOT_LOG2FC_GUIDE", inherits = TRUE)
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
) %>%
  mutate(title = gsub("_", "-", contrast))

panel_meta <- read_csv(PANEL_SOURCE, show_col_types = FALSE) %>%
  distinct(panel, PanelProtein, PanelGene, PanelUniProt, SomaKey)

read_contrast <- function(file, contrast, title) {
  df <- read_csv(file.path(DREAM_OUT_DIR, file), show_col_types = FALSE)
  n_col <- paste0("n_", unique(df$comparison_group)[1])

  df %>%
    mutate(
      contrast = contrast,
      contrast_title = title,
      n_comparison = .data[[n_col]]
    )
}

dream_res <- purrr::pmap_dfr(contrasts, read_contrast)

panel_res <- panel_meta %>%
  inner_join(dream_res, by = "SomaKey", relationship = "many-to-many")

if (nrow(panel_res) == 0L) {
  warning("No panel SomaKeys were present in the DREAM top tables; skipping panel volcano plots.")
  quit(save = "no", status = 0)
}

panel_res <- panel_res %>%
  group_by(panel, contrast) %>%
  mutate(
    panel_adj.P.Val = p.adjust(P.Value, method = "BH"),
    neg_log10_p = -log10(pmax(P.Value, .Machine$double.xmin)),
    fold_change = 2^logFC,
    effect_bin = case_when(
      P.Value < PLOT_P_THRESHOLD & logFC >= PLOT_LOG2FC_GUIDE  ~ "Higher than cutoff",
      P.Value < PLOT_P_THRESHOLD & logFC <= -PLOT_LOG2FC_GUIDE ~ "Lower than cutoff",
      TRUE                                                     ~ "Not highlighted"
    ),
    evidence = case_when(
      panel_adj.P.Val < 0.05 ~ "Panel FDR < 0.05",
      panel_adj.P.Val < 0.10 ~ "Panel FDR < 0.10",
      panel_adj.P.Val < 0.20 ~ "Panel FDR < 0.20",
      TRUE                   ~ "Other"
    ),
    evidence = factor(
      evidence,
      levels = c("Panel FDR < 0.05", "Panel FDR < 0.10", "Panel FDR < 0.20", "Other")
    ),
    effect_bin = factor(
      effect_bin,
      levels = c("Higher than cutoff", "Lower than cutoff", "Not highlighted")
    ),
    abs_log2fc_ge_10pct = abs(logFC) >= PLOT_LOG2FC_GUIDE,
    label = case_when(
      !is.na(PanelGene) & PanelGene != "" ~ PanelGene,
      !is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "" ~ EntrezGeneSymbol,
      !is.na(PanelProtein) & PanelProtein != "" ~ PanelProtein,
      TRUE ~ SomaKey
    )
  ) %>%
  ungroup()

dir.create(PANEL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(panel_res, file.path(PANEL_OUT_DIR, "dream_panel_subset_all_somakeys.csv"))

summary_tbl <- panel_res %>%
  group_by(panel, contrast, contrast_title, comparison_group, reference_group) %>%
  summarise(
    n_somakeys = n_distinct(SomaKey),
    n_panel_proteins = n_distinct(PanelProtein),
    n_comparison = first(n_comparison),
    n_HC = first(n_HC),
    group_mode = first(group_mode),
    min_baseline_age = first(min_baseline_age),
    model_formula = first(model_formula),
    panel_fdr_0_05 = sum(panel_adj.P.Val < 0.05, na.rm = TRUE),
    panel_fdr_0_10 = sum(panel_adj.P.Val < 0.10, na.rm = TRUE),
    panel_fdr_0_20 = sum(panel_adj.P.Val < 0.20, na.rm = TRUE),
    min_panel_fdr = min(panel_adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(PANEL_OUT_DIR, "dream_panel_fdr_summary.csv"))

panel_res %>%
  filter(panel_adj.P.Val < 0.10) %>%
  arrange(panel, contrast, panel_adj.P.Val) %>%
  select(panel, contrast, contrast_title, comparison_group, reference_group,
         PanelProtein, PanelGene, SomaKey, Target, EntrezGeneSymbol,
         logFC, fold_change, P.Value, panel_adj.P.Val, adj.P.Val,
         abs_log2fc_ge_10pct, n_comparison, n_HC, model_formula) %>%
  write_csv(file.path(PANEL_OUT_DIR, "dream_panel_hits_fdr10.csv"))

palette <- c(
  "Higher than cutoff" = plot_color("PLOT_COLOR_HIGHER", "#D55E00"),
  "Lower than cutoff" = plot_color("PLOT_COLOR_LOWER", "#0072B2"),
  "Not highlighted" = plot_color("PLOT_COLOR_NEUTRAL", "#a6adb7")
)

make_panel_plot <- function(df, plot_title = NULL) {
  signal_df <- df %>% filter(evidence != "Other")
  label_df <- df %>%
    group_by(panel, contrast) %>%
    arrange(P.Value, .by_group = TRUE) %>%
    mutate(label_me = panel_adj.P.Val < 0.20 | row_number() <= 6) %>%
    ungroup() %>%
    filter(label_me)

  p <- ggplot(df, aes(x = logFC, y = neg_log10_p)) +
    geom_hline(yintercept = -log10(PLOT_P_THRESHOLD), color = "grey78",
               linewidth = 0.35, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey78", linewidth = 0.35) +
    geom_vline(xintercept = c(-PLOT_LOG2FC_GUIDE, PLOT_LOG2FC_GUIDE),
               color = "grey65", linewidth = 0.3, linetype = "dotted") +
    geom_point(aes(color = effect_bin), alpha = 0.68, size = 1.55) +
    geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 3.0,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.35,
      seed = 1
    ) +
    labs(
      title = plot_title,
      x = "Model log2 fold change vs HC",
      y = "-log10 P value",
      color = NULL,
      shape = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = palette, breaks = names(palette),
                       drop = FALSE, guide = "none")

  if (nrow(signal_df) > 0L) {
    p <- p +
      geom_point(
        data = signal_df,
        aes(color = effect_bin, shape = evidence),
        alpha = 0.92, size = 2.35, stroke = 0.9
      ) +
      scale_shape_manual(
        values = c("Panel FDR < 0.05" = 16,
                   "Panel FDR < 0.10" = 17,
                   "Panel FDR < 0.20" = 1),
        drop = TRUE
      ) +
      guides(shape = guide_legend(
        override.aes = list(color = "black", alpha = 1, size = 2.7)
      ))
  }

  p
}

for (this_panel in unique(panel_res$panel)) {
  for (this_contrast in unique(panel_res$contrast)) {
    df <- panel_res %>%
      filter(panel == this_panel, contrast == this_contrast)
    if (nrow(df) == 0L) next
    title <- sprintf("%s panel | %s", tools::toTitleCase(this_panel), unique(df$contrast_title))
    plot_obj <- make_panel_plot(df, title)
    ggsave(
      file.path(PANEL_VOLCANO_PLOT_DIR, sprintf("volcano_%s_%s.png", this_panel, this_contrast)),
      plot_obj, width = 7.4, height = 5.6, dpi = 300
    )
    ggsave(
      file.path(PANEL_VOLCANO_PLOT_DIR, sprintf("volcano_%s_%s.pdf", this_panel, this_contrast)),
      plot_obj, width = 7.4, height = 5.6
    )
  }
}

message("Wrote panel subset tables to: ", PANEL_OUT_DIR)
message("Wrote panel volcano plots to: ", PANEL_VOLCANO_PLOT_DIR)
