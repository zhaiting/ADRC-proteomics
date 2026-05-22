#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(broom.mixed)
  library(broom)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
analysis_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  normalizePath(getwd())
}
if (!file.exists(file.path(analysis_dir, "correlation_model_setup.Rmd"))) {
  analysis_dir <- normalizePath("correlation")
}

root_dir <- normalizePath(file.path(analysis_dir, ".."))
out_base <- file.path(root_dir, "outputs")
out_dir <- file.path(analysis_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DISCOVERY_FDR <- 0.20
CORRELATION_FDR <- 0.20
DISCOVERY_TAG <- sprintf("fdr%02d", round(DISCOVERY_FDR * 100))

message("Loading corrected cohort/model objects from correlation_model_setup.Rmd ...")
old_wd <- getwd()
setwd(analysis_dir)
on.exit(setwd(old_wd), add = TRUE)
tmp_r <- tempfile(fileext = ".R")
invisible(knitr::purl("correlation_model_setup.Rmd", output = tmp_r, quiet = TRUE))
SKIP_RMD_PANEL_RUNS <- TRUE
source(tmp_r)

panel_file <- Sys.getenv("MANUSCRIPT_PANEL_SOMAKEYS", unset = "")
if (!nzchar(panel_file)) {
  stop("Set MANUSCRIPT_PANEL_SOMAKEYS to the local manuscript panel SomaKey CSV.")
}
panel_file <- normalizePath(panel_file, mustWork = FALSE)
if (!file.exists(panel_file)) {
  stop("Missing panel membership file: ", panel_file)
}

combined_layers <- tibble::tribble(
  ~discovery_source, ~base_dir,
  "All",             file.path(out_base, "dream_combined_baseAge55"),
  "Female",          file.path(out_base, "dream_combined_baseAge55_sexF"),
  "Male",            file.path(out_base, "dream_combined_baseAge55_sexM")
)

combined_contrasts <- tibble::tribble(
  ~contrast,             ~file,                              ~contrast_title,
  "AD_spectrum_vs_HC",   "toptable_AD_spectrum_vs_HC.csv",   "AD_spectrum vs HC",
  "LB_spectrum_vs_HC",   "toptable_LB_spectrum_vs_HC.csv",   "LB_spectrum vs HC",
  "PD_vs_HC",            "toptable_PD_vs_HC.csv",            "PD vs HC"
)

panel_meta <- readr::read_csv(panel_file, show_col_types = FALSE) %>%
  distinct(panel, PanelProtein, PanelGene, PanelUniProt, SomaKey)

read_layer_contrast <- function(discovery_source, base_dir,
                                contrast, file, contrast_title) {
  path <- file.path(base_dir, file)
  if (!file.exists(path)) {
    message("  [skip] missing combined toptable: ", path)
    return(tibble())
  }
  df <- readr::read_csv(path, show_col_types = FALSE)
  n_col <- paste0("n_", unique(df$comparison_group)[1])
  df %>%
    mutate(
      discovery_source = discovery_source,
      contrast = contrast,
      contrast_title = contrast_title,
      n_comparison = .data[[n_col]]
    )
}

combined_all_layers <- purrr::pmap_dfr(
  tidyr::crossing(combined_layers, combined_contrasts),
  read_layer_contrast
)

combined_panel_all <- panel_meta %>%
  inner_join(combined_all_layers, by = "SomaKey", relationship = "many-to-many") %>%
  group_by(panel, discovery_source, contrast) %>%
  mutate(
    panel_adj.P.Val = p.adjust(P.Value, method = "BH"),
    fold_change = 2^logFC,
    label = case_when(
      !is.na(PanelGene) & PanelGene != "" ~ PanelGene,
      !is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "" ~ EntrezGeneSymbol,
      !is.na(PanelProtein) & PanelProtein != "" ~ PanelProtein,
      TRUE ~ SomaKey
    )
  ) %>%
  ungroup()

sex_panel_dir <- file.path(out_base, "dream_combined_baseAge55", "panel_subsets_sexpanels")
dir.create(sex_panel_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(
  combined_panel_all,
  file.path(sex_panel_dir, "dream_panel_subset_all_somakeys.csv")
)

combined_panel_summary <- combined_panel_all %>%
  group_by(panel, discovery_source, contrast, contrast_title,
           comparison_group, reference_group) %>%
  summarise(
    n_somakeys = n_distinct(SomaKey),
    n_panel_proteins = n_distinct(PanelProtein),
    n_comparison = first(n_comparison),
    n_HC = first(n_HC),
    group_mode = first(group_mode),
    min_baseline_age = first(min_baseline_age),
    panel_fdr_threshold = sum(panel_adj.P.Val < DISCOVERY_FDR, na.rm = TRUE),
    min_panel_fdr = min(panel_adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(
  combined_panel_summary,
  file.path(sex_panel_dir, paste0("dream_panel_", DISCOVERY_TAG, "_summary.csv"))
)

biomarker_order <- c(
  "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
  "PTau217", "GFAP", "NfL"
)
correlation_panels <- c("All sex-adjusted", "Female-only", "Male-only")

collapse_unique <- function(x) {
  x <- sort(unique(x[!is.na(x) & x != ""]))
  if (length(x) == 0L) "-" else paste(x, collapse = " | ")
}

make_unique_display_labels <- function(labels, somakeys) {
  tibble(label = labels, SomaKey = somakeys) %>%
    mutate(soma_id = str_match(SomaKey, "^SL\\d+_([^_]+)_")[, 2]) %>%
    group_by(label) %>%
    mutate(display = if (n() > 1L) paste0(label, " (", soma_id, ")") else label) %>%
    ungroup() %>%
    pull(display)
}

stacked_column_order <- function(beta_mat) {
  if (ncol(beta_mat) <= 2L) return(colnames(beta_mat))
  corr_mat <- suppressWarnings(stats::cor(
    beta_mat,
    use = "pairwise.complete.obs",
    method = "pearson"
  ))
  corr_mat[!is.finite(corr_mat)] <- 0
  diag(corr_mat) <- 1
  hc <- stats::hclust(stats::as.dist(1 - corr_mat), method = "average")
  colnames(beta_mat)[hc$order]
}

make_candidates <- function(panel_name) {
  combined_panel_all %>%
    filter(panel == panel_name, !is.na(panel_adj.P.Val),
           panel_adj.P.Val < DISCOVERY_FDR) %>%
    group_by(panel, SomaKey) %>%
    summarise(
      label = first(label),
      EntrezGeneSymbol = first(EntrezGeneSymbol),
      min_q_disc = min(panel_adj.P.Val, na.rm = TRUE),
      discovery = paste(
        unique(paste0(discovery_source, ":", contrast_title)),
        collapse = " | "
      ),
      discovery_source = collapse_unique(discovery_source),
      selection_contrast = collapse_unique(contrast_title),
      .groups = "drop"
    ) %>%
    arrange(discovery_source, selection_contrast, min_q_disc, label)
}

fit_one_pair_sex_panel <- function(somakey,
                                   biom,
                                   df,
                                   panel_label,
                                   age_terms = AGE_TERMS,
                                   biomarker_scale = "log2") {
  missing_age <- setdiff(age_terms, names(df))
  if (length(missing_age) > 0L) {
    stop("Missing age adjustment column(s): ", paste(missing_age, collapse = ", "))
  }

  df_cc <- df %>%
    select(SubjectID, all_of(age_terms), Sex, group_combined,
           y = all_of(somakey), x = all_of(biom)) %>%
    filter(!is.na(y), !is.na(x), if_all(all_of(age_terms), ~ !is.na(.x)),
           !is.na(SubjectID), !is.na(group_combined))
  if (panel_label == "Female-only") {
    df_cc <- df_cc %>% filter(tolower(as.character(Sex)) %in% c("f", "female"))
  } else if (panel_label == "Male-only") {
    df_cc <- df_cc %>% filter(tolower(as.character(Sex)) %in% c("m", "male"))
  } else {
    df_cc <- df_cc %>% filter(!is.na(Sex))
  }
  if (nrow(df_cc) < 10L) return(NULL)

  if (biomarker_scale == "log2") {
    if (any(df_cc$x <= 0, na.rm = TRUE)) return(NULL)
    df_cc <- df_cc %>% mutate(x_model = log2(x))
  } else {
    df_cc <- df_cc %>% mutate(x_model = x)
  }
  if (n_distinct(df_cc$x_model) < 2L) return(NULL)

  df_cc$Sex <- droplevels(df_cc$Sex)
  df_cc$group_combined <- droplevels(df_cc$group_combined)

  fixed <- c("x_model", age_terms)
  if (panel_label == "All sex-adjusted" && nlevels(df_cc$Sex) >= 2) {
    fixed <- c(fixed, "Sex")
  }
  if (nlevels(df_cc$group_combined) >= 2) {
    fixed <- c(fixed, "group_combined")
  }

  has_rep <- any(table(df_cc$SubjectID) > 1)
  fml_txt <- paste("y ~", paste(fixed, collapse = " + "),
                   if (has_rep) "+ (1|SubjectID)" else "")
  fit <- try(if (has_rep) lmer(as.formula(fml_txt), data = df_cc, REML = FALSE)
             else         lm(as.formula(fml_txt), data = df_cc),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  tt <- if (inherits(fit, "lmerMod")) broom.mixed::tidy(fit, effects = "fixed")
        else                          broom::tidy(fit)
  tt %>%
    filter(term == "x_model") %>%
    transmute(
      SomaKey = somakey,
      biomarker = biom,
      correlation_panel = panel_label,
      biomarker_scale = biomarker_scale,
      beta = estimate,
      se = std.error,
      p = p.value,
      n_obs = nrow(df_cc),
      n_subj = n_distinct(df_cc$SubjectID),
      model = if (has_rep) "lmer" else "lm",
      model_formula = fml_txt
    )
}

source_colors <- c(
  "All" = "#4D4D4D",
  "Female" = "#CC79A7",
  "Male" = "#009E73",
  "Female | Male" = "#56B4E9",
  "All | Female" = "#999999",
  "All | Male" = "#999999",
  "All | Female | Male" = "#999999"
)
contrast_colors <- c(
  "AD_spectrum vs HC" = "#1B4F8C",
  "LB_spectrum vs HC" = "#8C2A00",
  "PD vs HC" = "#B8860B",
  "LB_spectrum vs HC | PD vs HC" = "#B0B0B0",
  "AD_spectrum vs HC | LB_spectrum vs HC" = "#B0B0B0",
  "AD_spectrum vs HC | PD vs HC" = "#B0B0B0",
  "AD_spectrum vs HC | LB_spectrum vs HC | PD vs HC" = "#B0B0B0"
)

plot_panel <- function(panel_name) {
  candidates <- make_candidates(panel_name)
  if (nrow(candidates) == 0L) {
    warning("No FDR<", DISCOVERY_FDR, " combined sex-panel candidates for ", panel_name)
    return(invisible(NULL))
  }

  readr::write_csv(
    candidates,
    file.path(analysis_dir,
              paste0("candidate_roster_", panel_name,
                     "_combined_sexpanels_", DISCOVERY_TAG, "_v2.csv"))
  )

  pair_grid <- tidyr::expand_grid(
    SomaKey = candidates$SomaKey,
    biomarker = biomarker_order,
    correlation_panel = correlation_panels
  )
  res_lmm <- purrr::pmap_dfr(
    pair_grid,
    ~ fit_one_pair_sex_panel(..1, ..2, samp_tbl_bio, ..3)
  )
  if (nrow(res_lmm) == 0L) stop("All sex-panel LMM fits returned NULL.")

  res_lmm <- res_lmm %>%
    group_by(correlation_panel) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    left_join(candidates %>%
                select(SomaKey, ProteinLabel = label, discovery,
                       min_q_disc, discovery_source, selection_contrast),
              by = "SomaKey") %>%
    mutate(
      ci_low = beta - 1.96 * se,
      ci_high = beta + 1.96 * se
    ) %>%
    select(ProteinLabel, SomaKey, discovery, min_q_disc,
           discovery_source, selection_contrast,
           correlation_panel, biomarker, biomarker_scale,
           beta, se, ci_low, ci_high, p, q, n_obs, n_subj,
           model, model_formula) %>%
    arrange(correlation_panel, q, desc(abs(beta)))

  out_tag <- paste(panel_name, "combined", "sexpanels", DISCOVERY_TAG, sep = "_")
  readr::write_csv(
    res_lmm,
    file.path(analysis_dir,
              paste0("protein_biomarker_associations_", out_tag, "_v2.csv"))
  )

  col_df <- candidates %>%
    mutate(
      display = make_unique_display_labels(label, SomaKey),
      display = factor(display, levels = display)
    )
  col_order <- as.character(col_df$display)

  mat_df <- res_lmm %>%
    left_join(col_df %>% select(SomaKey, display), by = "SomaKey") %>%
    mutate(
      row_id = paste(correlation_panel, biomarker, sep = "__"),
      row_id = factor(
        row_id,
        levels = as.vector(outer(correlation_panels, biomarker_order,
                                 paste, sep = "__"))
      ),
      display = factor(as.character(display), levels = col_order)
    )

  beta_mat <- mat_df %>%
    select(row_id, display, beta) %>%
    pivot_wider(names_from = display, values_from = beta) %>%
    arrange(row_id) %>%
    column_to_rownames("row_id") %>%
    as.matrix()

  q_mat <- mat_df %>%
    select(row_id, display, q) %>%
    pivot_wider(names_from = display, values_from = q) %>%
    arrange(row_id) %>%
    column_to_rownames("row_id") %>%
    as.matrix()

  column_order <- stacked_column_order(beta_mat)
  beta_mat <- beta_mat[, column_order, drop = FALSE]
  q_mat <- q_mat[, column_order, drop = FALSE]
  col_df <- col_df %>%
    mutate(display = factor(as.character(display), levels = column_order)) %>%
    arrange(display)

  tibble(
    display = column_order,
    column_order = seq_along(column_order),
    ordering_metric = "Pearson correlation distance on stacked beta matrix",
    linkage = "average"
  ) %>%
    left_join(
      col_df %>%
        mutate(display = as.character(display)) %>%
        select(display, SomaKey, discovery_source, selection_contrast, min_q_disc),
      by = "display"
    ) %>%
    readr::write_csv(
      file.path(analysis_dir,
                paste0("column_order_", out_tag, "_v2.csv"))
    )

  row_split <- factor(
    str_replace(rownames(beta_mat), "__.*$", ""),
    levels = correlation_panels
  )
  row_labels <- str_replace(rownames(beta_mat), "^.*__", "")
  star_mat <- ifelse(!is.na(q_mat) & q_mat < CORRELATION_FDR, "*", "")

  annotation_df <- col_df %>%
    select(display, discovery_source, selection_contrast) %>%
    column_to_rownames("display") %>%
    as.data.frame()

  source_levels <- unique(annotation_df$discovery_source)
  contrast_levels <- unique(annotation_df$selection_contrast)
  source_colors_use <- source_colors[source_levels]
  source_colors_use[is.na(source_colors_use)] <- "#999999"
  contrast_colors_use <- contrast_colors[contrast_levels]
  contrast_colors_use[is.na(contrast_colors_use)] <- "#B0B0B0"

  top_anno <- HeatmapAnnotation(
    `Discovery source` = annotation_df$discovery_source,
    `Selection contrast` = annotation_df$selection_contrast,
    col = list(
      `Discovery source` = source_colors_use,
      `Selection contrast` = contrast_colors_use
    ),
    show_legend = FALSE,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontface = "bold", fontsize = 7),
    simple_anno_size = unit(3, "mm")
  )

  mx <- max(0.30, abs(beta_mat), na.rm = TRUE)
  mx_display <- max(0.30, round(mx, 1))
  col_fun <- circlize::colorRamp2(c(-mx_display, 0, mx_display),
                                  c("#2166AC", "white", "#B2182B"))
  star_legend <- Legend(
    labels = sprintf("association FDR < %.2f", CORRELATION_FDR),
    title = NULL,
    type = "points",
    pch = "*",
    size = unit(2.2, "mm"),
    labels_gp = gpar(fontsize = 6.5)
  )
  beta_legend <- Legend(
    col_fun = col_fun,
    at = c(-mx_display, 0, mx_display),
    labels = c(sprintf("%.1f", -mx_display), "0", sprintf("%.1f", mx_display)),
    title = "beta",
    direction = "horizontal",
    title_position = "topcenter",
    legend_width = unit(18, "mm"),
    title_gp = gpar(fontface = "bold", fontsize = 6.5),
    labels_gp = gpar(fontsize = 6)
  )
  source_legend <- Legend(
    labels = source_levels,
    title = "Discovery source",
    legend_gp = gpar(fill = unname(source_colors_use), col = NA),
    grid_width = unit(2.2, "mm"),
    grid_height = unit(2.2, "mm"),
    nrow = 1,
    title_gp = gpar(fontface = "bold", fontsize = 6.5),
    labels_gp = gpar(fontsize = 6)
  )
  contrast_legend <- Legend(
    labels = contrast_levels,
    title = "Selection contrast",
    legend_gp = gpar(fill = unname(contrast_colors_use), col = NA),
    grid_width = unit(2.2, "mm"),
    grid_height = unit(2.2, "mm"),
    nrow = 1,
    title_gp = gpar(fontface = "bold", fontsize = 6.5),
    labels_gp = gpar(fontsize = 6)
  )
  bottom_legend <- packLegend(
    packLegend(beta_legend, star_legend, source_legend,
               direction = "horizontal", gap = unit(5, "mm")),
    contrast_legend,
    direction = "vertical",
    row_gap = unit(2, "mm")
  )
  panel_title <- recode(
    panel_name,
    bone = "Bone panel: plasma biomarker associations",
    senescence = "Senescence panel: plasma biomarker associations",
    .default = paste(panel_name, "panel: plasma biomarker associations")
  )

  ht <- Heatmap(
    beta_mat,
    name = "beta",
    col = col_fun,
    top_annotation = top_anno,
    width = unit(ncol(beta_mat) * 4.6, "mm"),
    height = unit(nrow(beta_mat) * 3.7, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    row_labels = row_labels,
    row_title_rot = 0,
    row_title_gp = gpar(fontface = "bold", fontsize = 8.5),
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 7.5),
    column_names_rot = 60,
    column_names_gp = gpar(fontsize = 7),
    column_title = panel_title,
    column_title_gp = gpar(fontface = "bold", fontsize = 11),
    rect_gp = gpar(col = "grey65", lwd = 0.4),
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(star_mat[i, j]) && star_mat[i, j] == "*") {
        grid.text("*", x, y, gp = gpar(fontsize = 8.5, fontface = "bold"))
      }
    }
  )

  plot_width <- if (ncol(beta_mat) > 8) 7.2 else 6.2
  plot_height <- 5.75
  out_stub <- paste("heatmap", panel_name, "combined", "sexpanels",
                    DISCOVERY_TAG, "v2", sep = "_")

  grDevices::pdf(file.path(out_dir, paste0(out_stub, ".pdf")),
                 width = plot_width, height = plot_height)
  draw(ht, heatmap_legend_side = "bottom",
       heatmap_legend_list = list(bottom_legend),
       show_annotation_legend = FALSE,
       legend_gap = unit(2, "mm"),
       padding = unit(c(3, 2, 3, 2), "mm"))
  grDevices::dev.off()

  grDevices::png(file.path(out_dir, paste0(out_stub, ".png")),
                 width = plot_width, height = plot_height,
                 units = "in", res = 300)
  draw(ht, heatmap_legend_side = "bottom",
       heatmap_legend_list = list(bottom_legend),
       show_annotation_legend = FALSE,
       legend_gap = unit(2, "mm"),
       padding = unit(c(3, 2, 3, 2), "mm"))
  grDevices::dev.off()

  message(sprintf(
    "%s: %d candidates; wrote %s.{pdf,png}",
    panel_name, nrow(candidates), file.path(out_dir, out_stub)
  ))
  invisible(list(candidates = candidates, associations = res_lmm))
}

bone_out <- plot_panel("bone")
senescence_out <- plot_panel("senescence")

message("Wrote FDR<", DISCOVERY_FDR, " combined sex-panel heatmaps to: ", out_dir)
