
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(patchwork)
})
detect_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
    return(normalizePath(file.path(dirname(script_path), "..")))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
  }
  normalizePath(getwd())
}
PROJECT_ROOT <- detect_project_root()

resolve_project_path <- function(path) {
  if (grepl("^/|^[A-Za-z]:", path)) {
    return(normalizePath(path, mustWork = FALSE))
  }
  normalizePath(file.path(PROJECT_ROOT, path), mustWork = FALSE)
}
getenv_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) default else value
}

RESULTS_IN_RAW <- getenv_or(
  "PATHWAY_RESULTS_DIRS",
  getenv_or("PATHWAY_RESULTS_IN",
            file.path(PROJECT_ROOT, "outputs", "dream_individual_baseAge55"))
)
RESULTS_DIRS <- strsplit(RESULTS_IN_RAW, ";", fixed = TRUE)[[1]]
RESULTS_DIRS <- trimws(RESULTS_DIRS[nzchar(trimws(RESULTS_DIRS))])
RESULTS_DIRS <- vapply(RESULTS_DIRS, resolve_project_path, character(1))
RESULTS_IN <- RESULTS_DIRS[[1]]
OUT_DIR <- resolve_project_path(getenv_or(
  "PATHWAY_OUT_DIR",
  file.path(PROJECT_ROOT, "pathway_analysis", "results")
))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "ranked_lists"), showWarnings = FALSE)
PLOT_DIR <- file.path(OUT_DIR, "plots")
EXPLORATORY_PLOT_DIR <- file.path(OUT_DIR, "exploratory_plots")
dir.create(PLOT_DIR,             showWarnings = FALSE)
dir.create(EXPLORATORY_PLOT_DIR, showWarnings = FALSE)

cat("Project root :", PROJECT_ROOT, "\n")
cat("Results dirs :\n")
cat(paste0("  - ", RESULTS_DIRS, collapse = "\n"), "\n")
cat("Output dir   :", OUT_DIR,      "\n\n")
DISEASE_ORDER <- c("AD_spectrum", "AD", "MCI", "LB_spectrum", "LBD", "PD", "MCI-PD")

comparison_disease <- function(x) {
  str_remove(x, "_vs_HC_(Combined|F|M)$")
}

comparison_stratum <- function(x) {
  case_when(
    str_detect(x, "_Combined$") ~ "All",
    str_detect(x, "_F$") ~ "F",
    str_detect(x, "_M$") ~ "M",
    TRUE ~ NA_character_
  )
}
GENE_AGG_METHOD <- getenv_or("PATHWAY_GENE_AGG", "median")
stopifnot(GENE_AGG_METHOD %in% c("median", "mean", "absmax"))
SYMBOL_EXCLUSIONS   <- c("None")
SYMBOL_EXCLUDE_PIPE <- TRUE

clean_symbol_mask <- function(symbols) {
  keep <- !is.na(symbols) & nzchar(symbols) & !(symbols %in% SYMBOL_EXCLUSIONS)
  if (SYMBOL_EXCLUDE_PIPE) keep <- keep & !grepl("\\|", symbols)
  keep
}

symbol_exclusion_audit <- tibble(source = character(),
                                 symbol = character(),
                                 reason = character(),
                                 n_rows = integer())

record_exclusions <- function(source_name, symbols) {
  if (length(symbols) == 0) return(invisible(NULL))
  reasons <- dplyr::case_when(
    is.na(symbols) | !nzchar(symbols)         ~ "na_or_empty",
    symbols %in% SYMBOL_EXCLUSIONS             ~ "excluded_literal",
    SYMBOL_EXCLUDE_PIPE & grepl("\\|", symbols) ~ "compound_pipe",
    TRUE ~ NA_character_
  )
  dropped <- !is.na(reasons)
  if (!any(dropped)) return(invisible(NULL))
  tibble(source = source_name,
         symbol = symbols[dropped],
         reason = reasons[dropped]) |>
    count(source, symbol, reason, name = "n_rows") -> add
  symbol_exclusion_audit <<- bind_rows(symbol_exclusion_audit, add)
}
cat("=== 1. Loading LMM results and building ranked lists ===\n")
cat(sprintf("  Gene aggregation: %s\n", GENE_AGG_METHOD))
cat(sprintf("  Symbol policy   : drop {%s}%s\n",
            paste(SYMBOL_EXCLUSIONS, collapse = ","),
            if (SYMBOL_EXCLUDE_PIPE) " + compound '|' symbols" else ""))

parse_filename <- function(fn) {
  base <- basename(fn)
  stem <- str_remove(base, "\\.csv$")
  stem <- str_remove(stem, "^toptable_")
  disease_raw <- str_remove(stem, "_vs_HC.*$")
  disease <- case_when(
    disease_raw == "MCI_PD" ~ "MCI-PD",
    TRUE ~ disease_raw
  )
  stratum <- if (str_detect(stem, "_Combined_")) "Combined"
             else if (str_detect(stem, "_Sex_F_")) "F"
             else if (str_detect(stem, "_Sex_M_")) "M"
             else NA_character_
  first <- tryCatch(read_csv(fn, show_col_types = FALSE, n_max = 1),
                    error = function(e) NULL)
  if (!is.null(first)) {
    if ("comparison_group" %in% names(first) &&
        !is.na(first$comparison_group[1]) &&
        nzchar(as.character(first$comparison_group[1]))) {
      disease <- as.character(first$comparison_group[1])
      disease <- ifelse(disease == "MCI_PD", "MCI-PD", disease)
    } else if ("contrast" %in% names(first) &&
               !is.na(first$contrast[1]) &&
               nzchar(as.character(first$contrast[1]))) {
      disease <- str_remove(as.character(first$contrast[1]), "_vs_HC$")
      disease <- ifelse(disease == "MCI_PD", "MCI-PD", disease)
    }
    if ("sex_stratum" %in% names(first) &&
        !is.na(first$sex_stratum[1]) &&
        nzchar(as.character(first$sex_stratum[1]))) {
      stratum <- case_when(
        tolower(as.character(first$sex_stratum[1])) %in% c("all", "combined") ~ "Combined",
        tolower(as.character(first$sex_stratum[1])) %in% c("female", "f") ~ "F",
        tolower(as.character(first$sex_stratum[1])) %in% c("male", "m") ~ "M",
        TRUE ~ stratum
      )
    }
  }
  if (is.na(stratum)) stratum <- "Combined"
  list(disease = disease, stratum = stratum,
       key = paste0(disease, "_vs_HC_", stratum))
}

protein_meta_path <- file.path(RESULTS_IN, "protein_metadata.csv")
protein_meta <- if (file.exists(protein_meta_path)) {
  read_csv(protein_meta_path, show_col_types = FALSE) |>
    select(Protein, EntrezGeneSymbol)
} else {
  NULL
}

pick_col <- function(df, candidates, default = NA) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(rep(default, nrow(df)))
  df[[hit[1]]]
}
per_comparison_files <- c(
  unlist(lapply(RESULTS_DIRS, function(results_dir) {
    c(
      list.files(results_dir, pattern = "_LMM_FDR_results\\.csv$", full.names = TRUE),
      list.files(file.path(results_dir, "per_comparison"),
                 pattern = "_LMM_FDR_results\\.csv$", full.names = TRUE),
      list.files(results_dir, pattern = "^toptable_.*_vs_HC\\.csv$", full.names = TRUE)
    )
  }), use.names = FALSE)
)
per_comparison_files <- unique(per_comparison_files)
if (length(per_comparison_files) == 0) {
  stop("No supported comparison files found in ", RESULTS_IN)
}
per_comparison_files <- map_dfr(per_comparison_files, function(fp) {
  m <- parse_filename(fp)
  tibble(file = fp, Disease = m$disease, Stratum = m$stratum)
}) |>
  mutate(Disease = factor(Disease, levels = DISEASE_ORDER),
         Stratum = factor(Stratum, levels = c("Combined", "F", "M"))) |>
  arrange(Disease, Stratum, file) |>
  pull(file)
cat("  found", length(per_comparison_files), "comparison files\n")

standardize_result_table <- function(fp) {
  raw <- read_csv(fp, show_col_types = FALSE)

  if ("z_Coefficient_Group" %in% names(raw)) {
    if (!"EntrezGeneSymbol" %in% names(raw)) {
      if (is.null(protein_meta)) {
        stop("Missing EntrezGeneSymbol and no protein_metadata.csv found for ", fp)
      }
      raw <- raw |> inner_join(protein_meta, by = "Protein")
    }
    return(tibble(
      EntrezGeneSymbol = as.character(raw$EntrezGeneSymbol),
      feature_id       = as.character(pick_col(raw, c("Protein", "SomaKey", "SomaId"))),
      Target           = as.character(pick_col(raw, c("Target", "TargetFullName"))),
      z                = as.numeric(raw$z_Coefficient_Group),
      effect           = as.numeric(pick_col(raw, c("Coefficient_Group", "logFC"))),
      p_value          = as.numeric(pick_col(raw, c("P_value_Group", "P.Value"))),
      padj             = as.numeric(pick_col(raw, c("P_value_Group_FDR_corrected", "adj.P.Val"))),
      source_file      = basename(fp)
    ))
  }

  if ("z.std" %in% names(raw)) {
    return(tibble(
      EntrezGeneSymbol = as.character(raw$EntrezGeneSymbol),
      feature_id       = as.character(pick_col(raw, c("SomaKey", "Protein", "SomaId"))),
      Target           = as.character(pick_col(raw, c("Target", "TargetFullName"))),
      z                = as.numeric(raw$z.std),
      effect           = as.numeric(pick_col(raw, c("logFC", "Coefficient_Group"))),
      p_value          = as.numeric(pick_col(raw, c("P.Value", "P_value_Group"))),
      padj             = as.numeric(pick_col(raw, c("adj.P.Val", "P_value_Group_FDR_corrected"))),
      source_file      = basename(fp)
    ))
  }

  stop("No supported z-statistic column found in ", fp)
}

aggregate_by_gene <- function(df, method) {
  agg <- df |>
    group_by(EntrezGeneSymbol) |>
    summarise(
      z          = switch(method,
                          median = median(z, na.rm = TRUE),
                          mean   = mean(z,   na.rm = TRUE),
                          absmax = z[which.max(abs(z))]),
      n_aptamers       = n(),
      all_feature_ids  = paste(feature_id, collapse = ";"),
      .groups          = "drop"
    )
  rep_row <- df |>
    group_by(EntrezGeneSymbol) |>
    slice_max(abs(z), n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(EntrezGeneSymbol,
           representative_feature_id = feature_id,
           representative_Target     = Target,
           representative_effect     = effect,
           representative_p_value    = p_value,
           representative_padj       = padj,
           source_file)
  agg |> left_join(rep_row, by = "EntrezGeneSymbol") |> arrange(desc(z))
}

ranked_lists <- list()
ranked_meta  <- tibble()
for (fp in per_comparison_files) {
  m   <- parse_filename(fp)
  raw_full <- standardize_result_table(fp) |> filter(!is.na(z))
  n_raw   <- nrow(raw_full)
  record_exclusions(paste0("ranked:", m$key), raw_full$EntrezGeneSymbol)
  keep_mask <- clean_symbol_mask(raw_full$EntrezGeneSymbol)
  raw       <- raw_full[keep_mask, ]
  n_excluded <- n_raw - nrow(raw)

  df <- aggregate_by_gene(raw, GENE_AGG_METHOD)

  ranks <- setNames(df$z, df$EntrezGeneSymbol)
  ranked_lists[[m$key]] <- ranks
  ranked_meta <- bind_rows(ranked_meta, tibble(
    Comparison              = m$key,
    Disease                 = m$disease,
    Stratum                 = m$stratum,
    n_raw_rows              = n_raw,
    n_symbol_excluded       = n_excluded,
    n_clean_rows            = nrow(raw),
    n_genes                 = length(ranks),
    n_aptamers_per_gene_med = median(df$n_aptamers),
    agg_method              = GENE_AGG_METHOD,
    file                    = basename(fp)))
  write_csv(df, file.path(OUT_DIR, "ranked_lists", paste0(m$key, "_ranked.csv")))
}
cat("\n  ranked-list summary:\n")
print(ranked_meta)
group_sizes_path <- file.path(RESULTS_IN, "group_sizes.csv")
group_sizes <- if (file.exists(group_sizes_path)) {
  read_csv(group_sizes_path, show_col_types = FALSE)
} else {
  NULL
}

extract_cohort_counts <- function(fp) {
  m <- parse_filename(fp)
  first <- read_csv(fp, show_col_types = FALSE, n_max = 1)

  if ("N_Subjects_Model" %in% names(first)) {
    return(tibble(Disease = m$disease, Stratum = m$stratum,
                  n_samples = first$N_Samples_Model[1],
                  n_subjects = first$N_Subjects_Model[1]))
  }

  disease_group <- str_replace(m$disease, "-", "_")
  local_group_sizes_path <- file.path(dirname(fp), "group_sizes.csv")
  local_group_sizes <- if (file.exists(local_group_sizes_path)) {
    read_csv(local_group_sizes_path, show_col_types = FALSE)
  } else {
    group_sizes
  }
  if (!is.null(local_group_sizes) &&
      all(c("group_var", "n_samples", "n_participants") %in% names(local_group_sizes))) {
    gs <- local_group_sizes |> filter(group_var %in% c(disease_group, "HC"))
    if (nrow(gs) > 0) {
      return(tibble(Disease = m$disease, Stratum = m$stratum,
                    n_samples = sum(gs$n_samples, na.rm = TRUE),
                    n_subjects = sum(gs$n_participants, na.rm = TRUE)))
    }
  }

  n_cols <- grep("^n_", names(first), value = TRUE)
  n_samples <- sum(unlist(first[n_cols]), na.rm = TRUE)
  tibble(Disease = m$disease, Stratum = m$stratum,
         n_samples = n_samples, n_subjects = n_samples)
}

cohort_n <- map_dfr(per_comparison_files, function(fp) {
  extract_cohort_counts(fp)
})

format_cohort_subtitle <- function(cohort_n, n_universe) {
  cohort_str <- cohort_n |>
    pivot_wider(id_cols = Disease,
                names_from = Stratum, values_from = n_subjects,
                names_glue = "n_{Stratum}") |>
    transmute(s = sprintf("%s %d/%d/%d",
                          Disease,
                          replace_na(get("n_Combined"), NA_integer_),
                          replace_na(get("n_F"), NA_integer_),
                          replace_na(get("n_M"), NA_integer_))) |>
    pull(s) |> paste(collapse = " | ")
  sprintf("Patients (All / F / M): %s    |    Gene universe: ~%s ranked proteins",
          cohort_str, formatC(n_universe, format = "d", big.mark = ","))
}
make_side_text_legend <- function(cohort_n, n_universe) {
  ord <- c("AD","MCI","LBD","PD","MCI-PD")
  rows <- cohort_n |>
    filter(Disease %in% ord) |>
    pivot_wider(id_cols = Disease, names_from = Stratum, values_from = n_subjects) |>
    arrange(factor(Disease, levels = ord)) |>
    transmute(s = sprintf("  %-7s %d / %d / %d",
                          Disease, Combined, `F`, M)) |> pull(s)
  txt <- paste0(
    "Patients (All / F / M)\n",
    paste(rows, collapse = "\n"),
    sprintf("\n\nGene universe\n  ~%s ranked proteins",
            formatC(n_universe, format = "d", big.mark = ",")),
    "\n\nSignificance\n  Black ring = padj < 0.1",
    "\n\nLayout\n  Diseases as facet columns,\n  strata (All/F/M) within"
  )
  ggplot() + theme_void() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
             family = "mono", size = 2.9, lineheight = 1.05, label = txt) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme(plot.margin = margin(8, 4, 4, 4))
}
cat("  cohort-size table:\n")
print(cohort_n)

cat("\n=== 2. Loading candidate gene-set sources ===\n")

read_gene_list <- function(path, label) {
  if (!nzchar(path)) {
    stop("Set ", label, " to a text file containing the corresponding gene list provided in the manuscript.")
  }
  path <- resolve_project_path(path)
  if (!file.exists(path)) stop("Gene list file not found: ", path)
  genes <- trimws(readLines(path, warn = FALSE))
  unique(genes[nzchar(genes)])
}

clean_geneset_members <- function(name, members) {
  record_exclusions(paste0("geneset:", name), members)
  members[clean_symbol_mask(members)]
}

gs_custom_raw <- list(
  `Self-defined_Bone` = read_gene_list(
    getenv_or("MANUSCRIPT_BONE_GENE_LIST", ""),
    "MANUSCRIPT_BONE_GENE_LIST"
  ),
  `Self-defined_Senescence` = read_gene_list(
    getenv_or("MANUSCRIPT_SENESCENCE_GENE_LIST", ""),
    "MANUSCRIPT_SENESCENCE_GENE_LIST"
  )
)
gs_custom <- Map(clean_geneset_members, names(gs_custom_raw), gs_custom_raw)
for (nm in names(gs_custom)) {
  n_drop <- length(gs_custom_raw[[nm]]) - length(gs_custom[[nm]])
  cat(sprintf("  %-18s  %d genes (%d excluded by symbol policy)\n",
              nm, length(gs_custom[[nm]]), n_drop))
}

hm_df <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
gs_hallmark <- split(hm_df$gene_symbol, hm_df$gs_name)
cat(sprintf("  Hallmark:           %d sets\n", length(gs_hallmark)))

re_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
gs_reactome_full <- split(re_df$gene_symbol, re_df$gs_name)
gs_reactome <- gs_reactome_full[lengths(gs_reactome_full) >= 15 &
                                 lengths(gs_reactome_full) <= 500]
cat(sprintf("  Reactome (15-500):  %d sets (of %d)\n",
            length(gs_reactome), length(gs_reactome_full)))

cgp_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP"),
  error = function(e) tibble()
)
saul_sen_mayo <- cgp_df |>
  filter(gs_name == "SAUL_SEN_MAYO") |>
  pull(gene_symbol) |>
  unique()
gs_published_sen <- if (length(saul_sen_mayo) > 0) {
  list(SAUL_SEN_MAYO = clean_geneset_members("SAUL_SEN_MAYO", saul_sen_mayo))
} else {
  list()
}
if (nrow(symbol_exclusion_audit) > 0) {
  write_csv(symbol_exclusion_audit |> arrange(source, reason, desc(n_rows)),
            file.path(OUT_DIR, "symbol_exclusion_audit.csv"))
  cat(sprintf("  symbol_exclusion_audit.csv : %d source-x-symbol rows (see file for breakdown)\n",
              nrow(symbol_exclusion_audit)))
}
cat("\n=== 3. Building Tier 1 focused library ===\n")
focused_hallmark_keep <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
reactome_bone_pat <- "BONE|OSSIFICATION|OSTEO|^REACTOME_COLLAGEN|^REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX|^REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION|^REACTOME_ECM_PROTEOGLYCANS|^REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS|^REACTOME_NCAM1_INTERACTIONS|^REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS|^REACTOME_SIGNALING_BY_BMP|^REACTOME_BMP_SIGNALING|^REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT|^REACTOME_SIGNALING_BY_WNT|^REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS|^REACTOME_SIGNALING_BY_TGF_BETA|^REACTOME_SIGNALING_BY_NOTCH$|^REACTOME_RUNX2|CHONDRO|SKELETAL"
reactome_senescence_pat <- "SENESCENCE|SASP|TELOMERE|HETEROCHROMATIN|^REACTOME_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS|^REACTOME_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS|^REACTOME_DNA_DAMAGE_BYPASS"
reactome_p53_pat <- "TP53|^REACTOME_REGULATION_OF_TP53|^REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53"

reactome_focused_names <- unique(c(
  grep(reactome_bone_pat,        names(gs_reactome_full), value = TRUE),
  grep(reactome_senescence_pat,  names(gs_reactome_full), value = TRUE),
  grep(reactome_p53_pat,         names(gs_reactome_full), value = TRUE)
))
gs_reactome_focused <- gs_reactome_full[intersect(reactome_focused_names,
                                                  names(gs_reactome_full))]
gs_reactome_focused <- gs_reactome_focused[
  lengths(gs_reactome_focused) >= 10 & lengths(gs_reactome_focused) <= 500]
gs_hallmark_focused <- gs_hallmark[intersect(focused_hallmark_keep, names(gs_hallmark))]
focused_pool <- c(gs_custom, gs_hallmark_focused, gs_published_sen, gs_reactome_focused)
focused_pool <- focused_pool[!duplicated(names(focused_pool))]
focused_source_map <- c(
  setNames(rep("Custom",    length(gs_custom)),          names(gs_custom)),
  setNames(rep("Hallmark",  length(gs_hallmark_focused)),names(gs_hallmark_focused)),
  setNames(rep("Published", length(gs_published_sen)),   names(gs_published_sen)),
  setNames(rep("Reactome",  length(gs_reactome_focused)),names(gs_reactome_focused)))
focused_source_map <- focused_source_map[names(focused_pool)]
theme_tag <- function(name, members) {
  if (name %in% c("Self-defined_Bone")) return("Bone-custom")
  if (name %in% c("Self-defined_Senescence")) return("Senescence-custom")
  if (name %in% c("SAUL_SEN_MAYO")) return("Senescence-published")
  if (str_detect(name, "BONE|OSSIFICATION|OSTEO|COLLAGEN|ECM|INTEGRIN|CHONDRO|SKELETAL|TGF_BETA|^HALLMARK_TGF|^HALLMARK_WNT|^HALLMARK_NOTCH|BMP|WNT")) return("Bone-biology")
  if (str_detect(name, "SENESCENCE|SASP|TELOMERE|HETEROCHROMATIN|TP53|P53|APOPTOSIS|DNA_REPAIR|DNA_DAMAGE|OXIDATIVE_STRESS")) return("Senescence-aging")
  if (str_detect(name, "INTERFERON|TNFA_SIGNALING|INFLAMMATORY|IL6_JAK|HYPOXIA|REACTIVE_OXYGEN|OXIDATIVE_PHOS|MYC_TARGETS")) return("Aging-context")
  if (str_detect(name, "ESTROGEN|ANDROGEN")) return("Sex-hormone")
  return("Other")
}
focused_themes <- mapply(theme_tag, names(focused_pool), focused_pool)

cat(sprintf("  Tier 1 focused library: %d pathways\n", length(focused_pool)))
cat(sprintf("    by source:  %s\n", paste(table(focused_source_map), names(table(focused_source_map)), sep="x", collapse=", ")))
cat(sprintf("    by theme:   %s\n", paste(table(focused_themes), names(table(focused_themes)), sep="x", collapse=", ")))
focused_index <- tibble(
  pathway = names(focused_pool),
  source  = focused_source_map[names(focused_pool)],
  theme   = focused_themes[names(focused_pool)],
  size    = lengths(focused_pool),
  members = sapply(focused_pool, paste, collapse = ";")
) |> arrange(theme, source, pathway)
write_csv(focused_index, file.path(OUT_DIR, "focused_pathway_library.csv"))
cat(sprintf("  saved focused library index -> %s\n",
            file.path(basename(OUT_DIR), "focused_pathway_library.csv")))
cat("\n=== 4. Building Tier 2 unified pool (supplementary) ===\n")
unified_pool <- c(gs_custom, gs_hallmark, gs_reactome)
unified_source_map <- c(
  setNames(rep("Custom",   length(gs_custom)),   names(gs_custom)),
  setNames(rep("Hallmark", length(gs_hallmark)), names(gs_hallmark)),
  setNames(rep("Reactome", length(gs_reactome)), names(gs_reactome)))
keep_idx <- lengths(unified_pool) >= 10 | names(unified_pool) %in% names(gs_custom)
unified_pool <- unified_pool[keep_idx]
unified_source_map <- unified_source_map[keep_idx]
cat(sprintf("  Tier 2 unified pool: %d pathways\n", length(unified_pool)))
cat("\n=== 5. Running fgsea ===\n")
set.seed(42)

run_fgsea_pool <- function(ranks, pool, source_map, comparison) {
  seed_text <- paste(comparison, length(pool), sep = "_")
  set.seed(42 + sum(utf8ToInt(seed_text)) %% 1000000)
  res <- fgsea::fgsea(pathways    = pool,
                      stats       = ranks,
                      minSize     = 10,
                      maxSize     = 500,
                      eps         = 0,
                      nPermSimple = 10000) |> as_tibble()
  if (nrow(res) == 0) return(tibble())
  res |>
    mutate(Comparison                  = comparison,
           Source                      = source_map[pathway],
           padj_within_contrast_tier   = p.adjust(pval, method = "BH"),
           padj                        = padj_within_contrast_tier,
           leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
    select(Comparison, Source, pathway, NES, ES, pval,
           padj_within_contrast_tier, padj, size, leadingEdge)
}
add_cross_contrast_fdr <- function(results_tbl) {
  if (nrow(results_tbl) == 0) return(results_tbl)
  results_tbl |>
    mutate(padj_across_all_contrasts = p.adjust(pval, method = "BH"))
}
cat("\n  Tier 1 (focused):\n")
focused_results_list <- list()
for (rname in names(ranked_lists)) {
  ranks <- ranked_lists[[rname]]
  res <- run_fgsea_pool(ranks, focused_pool, focused_source_map, rname)
  cat(sprintf("    %-25s -> %d pathways tested, %d sig (padj<0.1)\n",
              rname, nrow(res), sum(res$padj < 0.1, na.rm = TRUE)))
  focused_results_list[[rname]] <- res
}
focused_results <- bind_rows(focused_results_list) |>
  add_cross_contrast_fdr() |>
  mutate(Theme = focused_themes[pathway])
cat("\n  Tier 2 (unified, supplementary):\n")
unified_results_list <- list()
for (rname in names(ranked_lists)) {
  ranks <- ranked_lists[[rname]]
  res <- run_fgsea_pool(ranks, unified_pool, unified_source_map, rname)
  cat(sprintf("    %-25s -> %d pathways tested, %d sig (padj<0.1)\n",
              rname, nrow(res), sum(res$padj < 0.1, na.rm = TRUE)))
  unified_results_list[[rname]] <- res
}
unified_results <- bind_rows(unified_results_list) |>
  add_cross_contrast_fdr()
cat("\n=== 6. Writing output tables ===\n")

write_csv(focused_results, file.path(OUT_DIR, "gsea_focused.csv"))
focused_sig_within <- focused_results |>
  filter(padj_within_contrast_tier < 0.1) |>
  arrange(Comparison, padj_within_contrast_tier)
focused_sig_across <- focused_results |>
  filter(padj_across_all_contrasts < 0.1) |>
  arrange(Comparison, padj_across_all_contrasts)
write_csv(focused_sig_within,
          file.path(OUT_DIR, "gsea_focused_significant_within.csv"))
write_csv(focused_sig_across,
          file.path(OUT_DIR, "gsea_focused_significant_across_contrasts.csv"))
cat(sprintf("  gsea_focused.csv                                 -> %d rows (Tier 1 primary)\n",
            nrow(focused_results)))
cat(sprintf("  gsea_focused_significant_within.csv              -> %d rows (within-contrast/tier FDR<0.1)\n",
            nrow(focused_sig_within)))
cat(sprintf("  gsea_focused_significant_across_contrasts.csv    -> %d rows (across-contrasts FDR<0.1)\n",
            nrow(focused_sig_across)))

write_csv(unified_results, file.path(OUT_DIR, "gsea_unified.csv"))
unified_sig_within <- unified_results |>
  filter(padj_within_contrast_tier < 0.1) |>
  arrange(Comparison, padj_within_contrast_tier)
unified_sig_across <- unified_results |>
  filter(padj_across_all_contrasts < 0.1) |>
  arrange(Comparison, padj_across_all_contrasts)
write_csv(unified_sig_within,
          file.path(OUT_DIR, "gsea_unified_significant_within.csv"))
write_csv(unified_sig_across,
          file.path(OUT_DIR, "gsea_unified_significant_across_contrasts.csv"))
cat(sprintf("  gsea_unified.csv                                 -> %d rows (Tier 2 supplementary)\n",
            nrow(unified_results)))
cat(sprintf("  gsea_unified_significant_within.csv              -> %d rows (within-contrast/tier FDR<0.1)\n",
            nrow(unified_sig_within)))
cat(sprintf("  gsea_unified_significant_across_contrasts.csv    -> %d rows (across-contrasts FDR<0.1)\n",
            nrow(unified_sig_across)))
focused_sig <- focused_sig_within
unified_sig <- unified_sig_within
position_table <- bind_rows(
  focused_results |> group_by(Comparison) |>
    arrange(padj_within_contrast_tier, .by_group = TRUE) |>
    mutate(rank = row_number(), total = n(),
           percentile = 100 * (1 - (rank - 1) / total),
           Universe = "Tier1_focused") |>
    ungroup() |> filter(Source == "Custom"),
  unified_results |> group_by(Comparison) |>
    arrange(padj_within_contrast_tier, .by_group = TRUE) |>
    mutate(rank = row_number(), total = n(),
           percentile = 100 * (1 - (rank - 1) / total),
           Universe = "Tier2_unified") |>
    ungroup() |> filter(Source == "Custom")) |>
  mutate(padj = padj_within_contrast_tier) |>
  select(Universe, Comparison, pathway, NES, pval,
         padj_within_contrast_tier, padj_across_all_contrasts, padj,
         size, rank, total, percentile, leadingEdge) |>
  arrange(Comparison, pathway, Universe)
write_csv(position_table, file.path(OUT_DIR, "custom_panel_position.csv"))
cat(sprintf("  custom_panel_position.csv    -> %d rows (custom panel rank in both tiers)\n",
            nrow(position_table)))
cat("\n=== 7. Plotting ===\n")

theme_clean <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
plot_focused_dotplot <- function(df, themes, outfile) {
  d <- df |> mutate(Theme = themes[pathway],
                    label = str_replace_all(pathway, "_", " ") |> str_trunc(60))
  theme_order <- c("Bone-custom","Bone-biology",
                   "Senescence-custom","Senescence-published","Senescence-aging",
                   "Sex-hormone","Aging-context","Other")
  d$Theme <- factor(d$Theme, levels = theme_order)
  d <- d |> arrange(Theme, desc(abs(NES)))
  d$label <- factor(d$label, levels = unique(d$label))

  p <- ggplot(d, aes(x = Comparison, y = label,
                     fill = NES, size = -log10(pmax(padj, 1e-300)))) +
    geom_point(shape = 21, color = "grey20", stroke = 0.3) +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(2, 9), name = "-log10(padj)") +
    labs(title = "Tier 1 - focused bone / senescence / aging pathway analysis",
         subtitle = sprintf("%d pathways, FDR within focused universe", length(unique(d$pathway))),
         x = NULL, y = NULL) +
    theme_clean +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(angle = 0, face = "bold"),
          panel.spacing.y = unit(0.3, "lines"))
  ggsave(outfile, p,
         height = max(6, 0.32 * length(unique(d$label)) + 1.5),
         width = 10)
  cat("    wrote", basename(outfile), "\n")
}
plot_focused_dotplot(focused_results, focused_themes,
  file.path(EXPLORATORY_PLOT_DIR, "tier1_focused_dotplot.pdf"))
plot_focused_heatmap <- function(df, themes, outfile) {
  d <- df
  mat <- d |> select(Comparison, pathway, NES) |>
    pivot_wider(names_from = Comparison, values_from = NES, values_fill = 0) |>
    column_to_rownames("pathway") |> as.matrix()
  pad <- d |> select(Comparison, pathway, padj) |>
    pivot_wider(names_from = Comparison, values_from = padj, values_fill = 1) |>
    column_to_rownames("pathway") |> as.matrix()
  pad <- pad[rownames(mat), colnames(mat), drop = FALSE]
  src_v <- themes[rownames(mat)]

  rownames(mat) <- str_replace_all(rownames(mat), "_", " ") |> str_trunc(60)
  rownames(pad) <- rownames(mat)
  theme_order <- c("Bone-custom","Bone-biology",
                   "Senescence-custom","Senescence-published","Senescence-aging",
                   "Sex-hormone","Aging-context","Other")
  ord <- order(factor(src_v, levels = theme_order),
               -apply(abs(mat), 1, max))
  mat <- mat[ord, , drop = FALSE]; pad <- pad[ord, , drop = FALSE]; src_v <- src_v[ord]

  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {
    col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5),
                                    c("#2b6cb0", "white", "#c1432b"))
    cell_fun <- function(j, i, x, y, w, h, fill) {
      if (!is.na(pad[i, j]) && pad[i, j] < 0.1) {
        grid::grid.text("*", x, y,
                        gp = grid::gpar(fontsize = 12, col = "black"))
      }
    }
    theme_palette <- c(`Bone-custom`="#d97706",`Bone-biology`="#f59e0b",
                       `Senescence-custom`="#7c3aed",`Senescence-published`="#a855f7",
                       `Senescence-aging`="#c084fc",
                       `Sex-hormone`="#ec4899",`Aging-context`="#0ea5e9",`Other`="grey70")
    row_anno <- ComplexHeatmap::rowAnnotation(
      Theme = src_v,
      col = list(Theme = theme_palette))
    ht <- ComplexHeatmap::Heatmap(
      mat, name = "NES", col = col_fun, cell_fun = cell_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      left_annotation = row_anno,
      row_names_max_width = grid::unit(20, "cm"),
      column_title = "Tier 1 - focused pathway NES across strata",
      column_title_gp = grid::gpar(fontface = "bold"))
    pdf(outfile, height = max(5, 0.27 * nrow(mat) + 2), width = 10)
    ComplexHeatmap::draw(ht)
    dev.off()
    cat("    wrote", basename(outfile), "\n")
  }
}
plot_focused_heatmap(focused_results, focused_themes,
  file.path(EXPLORATORY_PLOT_DIR, "tier1_focused_heatmap.pdf"))
plot_custom_position <- function(df_position, outfile) {
  d <- df_position |> filter(Universe == "Tier1_focused") |>
    mutate(label = sprintf("%s [rank %d/%d, %.1f%%ile]",
                           pathway, rank, total, percentile))
  if (nrow(d) == 0) return(invisible())
  p <- ggplot(d, aes(x = NES, y = label, color = padj < 0.1)) +
    geom_segment(aes(x = 0, xend = NES, yend = label), linewidth = 0.6) +
    geom_point(aes(size = -log10(pmax(pval, 1e-300))), alpha = 0.9) +
    facet_wrap(~ Comparison, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c(`TRUE` = "#c1432b", `FALSE` = "grey50"),
                       name = "padj < 0.1") +
    scale_size_continuous(range = c(2, 7), name = "-log10(pval)") +
    labs(title = "Custom panel rank within Tier 1 focused universe",
         x = "NES", y = NULL) +
    theme_clean
  ggsave(outfile, p, height = 5, width = 10)
  cat("    wrote", basename(outfile), "\n")
}
plot_custom_position(position_table,
  file.path(EXPLORATORY_PLOT_DIR, "tier1_custom_panel_position.pdf"))
plot_unified_top_per_stratum <- function(df, top_n = 15, outfile) {
  d <- df |>
    group_by(Comparison) |>
    slice_min(padj, n = top_n) |>
    ungroup() |>
    mutate(label = str_replace_all(pathway, "_", " ") |> str_trunc(70))
  if (nrow(d) == 0) return(invisible())
  p <- ggplot(d, aes(x = NES, y = reorder(label, NES), fill = Source)) +
    geom_col() +
    facet_wrap(~ Comparison, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(Custom="#d97706", Hallmark="#7c3aed", Reactome="#2563eb")) +
    labs(title = paste("Tier 2 (supplementary) - top", top_n,
                       "per stratum, unified Hallmark + Reactome universe"),
         x = "NES", y = NULL) +
    theme_clean
  ggsave(outfile, p, height = 9, width = 11)
  cat("    wrote", basename(outfile), "\n")
}
plot_unified_top_per_stratum(unified_results, top_n = 15,
  file.path(EXPLORATORY_PLOT_DIR, "tier2_unified_top15.pdf"))
plot_unified_heatmap <- function(df, top_n_per_stratum = 20, outfile) {
  top_pathways <- df |>
    group_by(Comparison) |>
    slice_min(padj, n = top_n_per_stratum) |>
    pull(pathway) |> unique()
  d <- df |> filter(pathway %in% top_pathways)
  if (nrow(d) == 0) return(invisible())
  mat <- d |> select(Comparison, pathway, NES) |>
    pivot_wider(names_from = Comparison, values_from = NES, values_fill = 0) |>
    column_to_rownames("pathway") |> as.matrix()
  pad <- d |> select(Comparison, pathway, padj) |>
    pivot_wider(names_from = Comparison, values_from = padj, values_fill = 1) |>
    column_to_rownames("pathway") |> as.matrix()
  pad <- pad[rownames(mat), colnames(mat), drop = FALSE]
  src_v <- d |> distinct(pathway, Source) |> deframe()
  src_v <- src_v[rownames(mat)]
  rownames(mat) <- str_replace_all(rownames(mat), "_", " ") |> str_trunc(60)
  rownames(pad) <- rownames(mat)
  ord <- order(-rowSums(abs(mat)))
  mat <- mat[ord, , drop = FALSE]; pad <- pad[ord, , drop = FALSE]; src_v <- src_v[ord]

  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {
    col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5),
                                    c("#2b6cb0", "white", "#c1432b"))
    cell_fun <- function(j, i, x, y, w, h, fill) {
      if (!is.na(pad[i, j]) && pad[i, j] < 0.1) {
        grid::grid.text("*", x, y,
                        gp = grid::gpar(fontsize = 12, col = "black"))
      }
    }
    row_anno <- ComplexHeatmap::rowAnnotation(
      Source = src_v,
      col = list(Source = c(Custom="#d97706", Hallmark="#7c3aed", Reactome="#2563eb")))
    ht <- ComplexHeatmap::Heatmap(
      mat, name = "NES", col = col_fun, cell_fun = cell_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      left_annotation = row_anno,
      row_names_max_width = grid::unit(20, "cm"),
      column_title = "Tier 2 (supplementary) - unified-pool NES across strata",
      column_title_gp = grid::gpar(fontface = "bold"))
    pdf(outfile, height = max(5, 0.27 * nrow(mat) + 2), width = 10)
    ComplexHeatmap::draw(ht)
    dev.off()
    cat("    wrote", basename(outfile), "\n")
  }
}
plot_unified_heatmap(unified_results, top_n_per_stratum = 20,
  file.path(EXPLORATORY_PLOT_DIR, "tier2_unified_heatmap.pdf"))
interpretation_theme_patterns <- list(
  Vesicle_trafficking   = "VESICLE|TRAFFICK|CLATHRIN|ENDOCYTOSIS|^REACTOME_RAB|MEMBRANE_TRAFFIC|GOLGI|COPII|ENDOSOM|LYSOSOM",
  Calcium_signaling     = "CALCIUM|^REACTOME_IP3|FCERI|CA2_|RAS_ACTIVATION_UPON_CA2|CA_2_INFLUX|CASR",
  Ubiquitin_proteasome  = "UBIQUI|PROTEASOM|DEUBIQ",
  Interferon_antiviral  = "INTERFERON|ANTIVIRAL|^REACTOME_IFN|^REACTOME_ZBP1|HALLMARK_INTERFERON",
  TP53_apoptosis        = "TP53|^HALLMARK_P53|HALLMARK_APOPTOSIS|^REACTOME_REGULATION_OF_TP53|^REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53",
  Mitochondrial         = "MITOCH|OXIDATIVE_PHOSPHO|RESPIRATORY_ELECTRON|TCA_CYCLE|CITRIC_ACID|FATTY_ACID_BETA",
  Inflammation          = "INFLAMMATORY|TNFA_SIGNALING|^REACTOME_TNF|IL6_JAK|HALLMARK_HYPOXIA",
  PI3K_AKT              = "PI3K|AKT_MTOR",
  Translation           = "TRANSLATION|RIBOSOM|TRNA",
  Bone_ECM              = "BONE|OSSIFICATION|OSTEO|COLLAGEN|ECM|INTEGRIN|RUNX2|BMP|^REACTOME_SIGNALING_BY_TGF|^REACTOME_SIGNALING_BY_WNT|^REACTOME_SIGNALING_BY_NOTCH",
  Senescence            = "SENESCENCE|SASP|TELOMERE|HETEROCHROMATIN|SAUL_SEN_MAYO"
)
tag_interpretation_theme <- function(name) {
  for (theme in names(interpretation_theme_patterns)) {
    if (grepl(interpretation_theme_patterns[[theme]], name, ignore.case = FALSE,
              perl = TRUE)) return(theme)
  }
  return("Other")
}
unified_results_themed <- unified_results |>
  mutate(Interpretation_Theme = sapply(pathway, tag_interpretation_theme))
write_csv(unified_results_themed, file.path(OUT_DIR, "gsea_unified_themed.csv"))
cat(sprintf("\n  gsea_unified_themed.csv      -> %d rows (Tier 2 with interpretation-theme tag)\n",
            nrow(unified_results_themed)))

unified_sig_themed <- unified_results_themed |> filter(padj < 0.1) |>
  arrange(Interpretation_Theme, Comparison, padj)
write_csv(unified_sig_themed,
          file.path(OUT_DIR, "gsea_unified_significant_themed.csv"))
cat(sprintf("  gsea_unified_significant_themed.csv -> %d rows (themed sig hits)\n",
            nrow(unified_sig_themed)))
plot_interpretation_themes_dotplot <- function(df_themed, outfile) {
  d <- df_themed |> filter(Interpretation_Theme != "Other", padj < 0.25) |>
    mutate(label = str_replace_all(pathway, "_", " ") |> str_trunc(60),
           Interpretation_Theme = factor(Interpretation_Theme,
             levels = c("Bone_ECM","Senescence","TP53_apoptosis",
                        "Vesicle_trafficking","Calcium_signaling",
                        "Ubiquitin_proteasome","Interferon_antiviral",
                        "Mitochondrial","Inflammation","Translation","PI3K_AKT")))
  if (nrow(d) == 0) return(invisible())
  d <- d |> arrange(Interpretation_Theme, desc(abs(NES)))
  d$label <- factor(d$label, levels = unique(d$label))
  p <- ggplot(d, aes(x = Comparison, y = label,
                     fill = NES, size = -log10(pmax(padj, 1e-300)))) +
    geom_point(shape = 21, color = "grey20", stroke = 0.3) +
    facet_grid(Interpretation_Theme ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(2, 9), name = "-log10(padj)") +
    labs(title = "Tier 2 (supplementary) - pathways grouped by interpretation theme",
         subtitle = "Pathways shown if padj<0.25 in any stratum; from unified Hallmark+Reactome pool",
         x = NULL, y = NULL) +
    theme_clean +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(angle = 0, face = "bold"),
          panel.spacing.y = unit(0.3, "lines"))
  ggsave(outfile, p,
         height = max(8, 0.28 * length(unique(d$label)) + 2),
         width = 11)
  cat("    wrote", basename(outfile), "\n")
}
plot_interpretation_themes_dotplot(unified_results_themed,
  file.path(EXPLORATORY_PLOT_DIR, "tier2_interpretation_themes_dotplot.pdf"))

cat("\n=== 7.5  Decomposing Self-defined_Senescence by hallmark sub-category ===\n")

classify_sen_hallmark <- function(g) {
  if (g %in% c("CDKN1A","CDKN2A","CDKN2B","CDKN1B","CDKN1C","CDKN2C","CDKN2D",
               "TP53","TP63","TP73","RB1","E2F1","CCNA2","CCNB1","CCND1","CCND2",
               "CCNE2","CCNG1","MKI67","CDK4","CDK6","BTG2","ETS2","EDA2R"))
    return("H1 Cell-cycle / arrest")
  if (g %in% c("H2AFX","H2AFJ","TP53BP1","MDM2","GADD45A","PHLDA3","PMAIP1","PPM1D",
               "DDB2","DDIT4","SESN2","TP53INP1","ZMAT3","FDXR","TREX2"))
    return("H2 DNA damage / DDR")
  if (g %in% c("LMNB1","LMNB2","HMGB1","HMGB2"))
    return("H3 Nuclear / chromatin")
  if (g %in% c("BCL2L1","FAS","TNFRSF1A","TNFRSF1B","TNFRSF10B","TNFRSF11B"))
    return("H4 Anti-apoptotic / death receptors")
  if (str_detect(g, "^IL[0-9]") ||
      g %in% c("TNF","TNFAIP2","TNFAIP6","CSF1","CSF2","MIF","IFNG","IFI44","PTGS2"))
    return("H5a SASP cytokines")
  if (str_detect(g, "^CXCL|^CCL[0-9]") || g %in% c("RARRES2"))
    return("H5b SASP chemokines")
  if (str_detect(g, "^MMP[0-9]|^TIMP[0-9]") ||
      g %in% c("ADAM19","HPSE"))
    return("H5c SASP MMPs / TIMPs")
  if (str_detect(g, "^IGFBP[0-9]|^FGF[0-9]") ||
      g %in% c("HGF","TGFB1","TGFA","NGF","PGF","PIGF","VEGFA","VEGFB","VEGFC","VEGFD",
                "EGF","EGFR","AREG","EREG","BMP2","BMP6","KITLG","ANG","ANGPTL4",
                "PDGFA","PTN","EDN1"))
    return("H5d SASP growth factors / IGFBPs")
  if (g %in% c("SERPINE1","SERPINE2","SERPINB2","SERPINB7","SERPINB4",
                "PLAU","PLAUR","PLAT","PLG"))
    return("H5e SASP serpins / fibrinolytic")
  if (g %in% c("THBS1","GSN","CHI3L1","GDF15","SPP1","PAPPA","ESM1","ANXA1"))
    return("H5f SASP secreted misc")
  if (g %in% c("GLB1","CTSB","CTSL"))
    return("H6 Lysosomal / SA-betagal")
  if (g %in% c("NAMPT","ALDH1A3","PRKN","NOS3","SELENOH","SIRT1"))
    return("H7 Metabolic / mitochondrial / ROS")
  if (g %in% c("DPP4","CD82","CD207","AGER","EPHA2","PECAM1","CD163"))
    return("H8 Cell-surface markers")
  if (str_detect(g, "^COL[0-9]") ||
      g %in% c("FN1","DCN","LOX","FBN1","SPARC","POSTN","HAS2","PCOLCE","LRRC17",
                "MFAP5","EFEMP1","FSTL1","CCDC80","ANKRD1","CRIM1","EDIL3","GPC1"))
    return("ECM / matrix")
  if (str_detect(g, "^KRT[0-9]"))
    return("Non-can: Keratin tissue markers")
  if (str_detect(g, "^CDH[0-9]"))
    return("Non-can: Cadherin tissue markers")
  if (str_detect(g, "^KLK[0-9]") ||
      g %in% c("CASP14","DNASE1L2","LAD1","LGALS7","CARD18","SPRR4","ALOX12B"))
    return("Non-can: Skin / keratinization")
  if (g %in% c("MX1","MX2","ISG15","IFIH1","IRF1","IRF7","GBP2","STAT1","RTP4",
                "APOBEC3G","TLR4","ZC3HAV1","USP18","SLC25A28"))
    return("Non-can: Antiviral / interferon")
  if (g %in% c("VWF","TEK","CDH5","FLT1","LYVE1"))
    return("Non-can: Vascular / endothelial")
  if (g %in% c("WIF1","WNT5A","WNT5B","WNT7B","SFRP1","CTNNB1","ROR1"))
    return("Non-can: WNT signaling")
  if (g %in% c("HSPA8","DNAJB1","HSPB6","CANX","UCHL1"))
    return("Non-can: Proteostasis / chaperone")
  if (g %in% c("APP","APOE"))
    return("Non-can: AD-specific")
  if (g %in% c("PSG4","PSG5"))
    return("Non-can: Pregnancy-specific")
  if (g %in% c("MAPK1","MAPK3","MAPK14","BRAF","RRAS"))
    return("Non-can: MAPK / stress signaling")
  if (g %in% c("ICAM1","VCAM1","SELPLG","SIRPA","TYROBP","ITGA2","ITGAL","CD44","CD47",
                "CD55","CD9","CD68","CD8A","CD8B","CD4","CD3D","MS4A1","FOXP3",
                "HLA-E","HLA-DQB1","HLA-DRB1","VSIR","FOS"))
    return("Non-can: Immune surface / adhesion")
  return("Other (uncategorized)")
}
sen_custom <- gs_custom[["Self-defined_Senescence"]]
sen_categories <- tibble(gene = sen_custom,
                         category = sapply(sen_custom, classify_sen_hallmark))

decomp_rows <- list()
for (rname in names(ranked_lists)) {
  ranks <- ranked_lists[[rname]]
  d <- tibble(gene = names(ranks), z = as.numeric(ranks)) |>
    inner_join(sen_categories, by = "gene") |>
    mutate(Comparison = rname)
  decomp_rows[[rname]] <- d
}
decomp_df <- bind_rows(decomp_rows)

decomp_summary <- decomp_df |>
  group_by(Comparison, category) |>
  summarise(n_in_assay = n(),
            mean_z     = mean(z),
            median_z   = median(z),
            n_up       = sum(z > 0),
            n_down     = sum(z < 0),
            pct_up     = 100 * sum(z > 0) / n(),
            .groups = "drop") |>
  arrange(category, Comparison)

write_csv(decomp_summary,
          file.path(OUT_DIR, "senescence_decomposition.csv"))
cat(sprintf("  wrote senescence_decomposition.csv (%d rows)\n", nrow(decomp_summary)))
plot_sen_decomp <- function(df, outfile) {
  cat_order <- c(
    "H1 Cell-cycle / arrest", "H2 DNA damage / DDR", "H3 Nuclear / chromatin",
    "H4 Anti-apoptotic / death receptors",
    "H5a SASP cytokines", "H5b SASP chemokines", "H5c SASP MMPs / TIMPs",
    "H5d SASP growth factors / IGFBPs", "H5e SASP serpins / fibrinolytic",
    "H5f SASP secreted misc",
    "H6 Lysosomal / SA-betagal", "H7 Metabolic / mitochondrial / ROS",
    "H8 Cell-surface markers",
    "ECM / matrix",
    "Non-can: Antiviral / interferon", "Non-can: Proteostasis / chaperone",
    "Non-can: WNT signaling", "Non-can: MAPK / stress signaling",
    "Non-can: Vascular / endothelial", "Non-can: Immune surface / adhesion",
    "Non-can: AD-specific", "Non-can: Skin / keratinization",
    "Non-can: Keratin tissue markers", "Non-can: Cadherin tissue markers",
    "Non-can: Pregnancy-specific",
    "Other (uncategorized)")
  d <- df |> mutate(category = factor(category, levels = rev(cat_order)),
                    label = sprintf("%s (n=%d)", category, n_in_assay))
  p <- ggplot(d, aes(x = Comparison, y = category, fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.2f\n(n=%d)", mean_z, n_in_assay)),
              size = 3, color = "black") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-0.6, 0.6),
                         oob = scales::squish, name = "mean z") +
    labs(title = "Decomposition of Self-defined senescence panel by hallmark sub-category",
         subtitle = "Mean z-statistic across genes in each sub-category. Cells show mean z and n genes measured.",
         x = NULL, y = NULL) +
    theme_clean +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(face = "bold"))
  ggsave(outfile, p, height = 7, width = 9)
  cat("    wrote", basename(outfile), "\n")
}
plot_sen_decomp(decomp_summary,
  file.path(EXPLORATORY_PLOT_DIR, "senescence_decomposition.pdf"))
plot_sen_decomp_bars <- function(df, outfile) {
  cat_order <- c(
    "H1 Cell-cycle / arrest", "H2 DNA damage / DDR", "H3 Nuclear / chromatin",
    "H4 Anti-apoptotic / death receptors",
    "H5a SASP cytokines", "H5b SASP chemokines", "H5c SASP MMPs / TIMPs",
    "H5d SASP growth factors / IGFBPs", "H5e SASP serpins / fibrinolytic",
    "H5f SASP secreted misc",
    "H6 Lysosomal / SA-betagal", "H7 Metabolic / mitochondrial / ROS",
    "H8 Cell-surface markers",
    "ECM / matrix",
    "Non-can: Antiviral / interferon", "Non-can: Proteostasis / chaperone",
    "Non-can: WNT signaling", "Non-can: MAPK / stress signaling",
    "Non-can: Vascular / endothelial", "Non-can: Immune surface / adhesion",
    "Non-can: AD-specific", "Non-can: Skin / keratinization",
    "Non-can: Keratin tissue markers", "Non-can: Cadherin tissue markers",
    "Non-can: Pregnancy-specific",
    "Other (uncategorized)")
  d <- df |> mutate(category = factor(category, levels = rev(cat_order)),
                    pct_down = 100 - pct_up)
  d_long <- d |>
    select(Comparison, category, n_in_assay, pct_up, pct_down) |>
    pivot_longer(c(pct_up, pct_down), names_to = "direction", values_to = "pct")
  p <- ggplot(d_long, aes(x = pct, y = category, fill = direction)) +
    geom_col(width = 0.8) +
    facet_wrap(~ Comparison, ncol = 3) +
    scale_fill_manual(values = c(pct_up = "#c1432b", pct_down = "#2b6cb0"),
                      labels = c(pct_up = "z > 0 (UP in disease)",
                                 pct_down = "z < 0 (DOWN in disease)"),
                      name = NULL) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "grey40") +
    labs(title = "Self-defined senescence panel decomposition - % genes UP vs DOWN per sub-category",
         x = "% of genes in sub-category", y = NULL) +
    theme_clean +
    theme(axis.text.y = element_text(face = "bold"),
          legend.position = "top")
  ggsave(outfile, p, height = 7, width = 12)
  cat("    wrote", basename(outfile), "\n")
}
plot_sen_decomp_bars(decomp_summary,
  file.path(EXPLORATORY_PLOT_DIR, "senescence_decomposition_bars.pdf"))
plot_sen_decomp_grid <- function(df, cohort_n, n_universe, outfile) {
  cat_order_full <- c(
    "H1 Cell-cycle / arrest", "H2 DNA damage / DDR", "H3 Nuclear / chromatin",
    "H4 Anti-apoptotic / death receptors",
    "H5a SASP cytokines", "H5b SASP chemokines", "H5c SASP MMPs / TIMPs",
    "H5d SASP growth factors / IGFBPs", "H5e SASP serpins / fibrinolytic",
    "H5f SASP secreted misc",
    "H6 Lysosomal / SA-betagal", "H7 Metabolic / mitochondrial / ROS",
    "H8 Cell-surface markers", "ECM / matrix",
    "Non-can: Antiviral / interferon", "Non-can: Proteostasis / chaperone",
    "Non-can: WNT signaling", "Non-can: MAPK / stress signaling",
    "Non-can: Vascular / endothelial", "Non-can: Immune surface / adhesion",
    "Non-can: AD-specific", "Non-can: Skin / keratinization",
    "Non-can: Keratin tissue markers", "Non-can: Cadherin tissue markers",
    "Non-can: Pregnancy-specific", "Other (uncategorized)")
  pretty <- function(x) {
    x |> str_replace("\\u03B2", "beta") |>
      str_replace("^H[0-9][a-z]? ", "") |>
      str_replace("^Non-can: ", "") |>
      str_replace("^B[0-9]+ ", "")
  }
  d <- df |>
    mutate(category_full = factor(category, levels = cat_order_full),
           category_pretty = factor(pretty(as.character(category_full)),
                                    levels = unique(pretty(cat_order_full))),
           Disease = comparison_disease(Comparison),
           Stratum = comparison_stratum(Comparison),
           Disease = factor(Disease, levels = DISEASE_ORDER),
           Stratum = factor(Stratum, levels = c("All","F","M")))

  p_main <- ggplot(d, aes(x = Stratum, y = forcats::fct_rev(category_pretty),
                          fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%+.2f", mean_z)),
              color = "grey15", size = 2.4) +
    facet_grid(. ~ Disease, switch = "y") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-1, 1),
                         oob = scales::squish, name = "mean z") +
    labs(title = "Self-defined senescence panel - gene-level decomposition by SenNet hallmark",
         x = NULL, y = NULL) +
    theme_clean +
    theme(strip.text.x = element_text(face = "bold"),
          axis.text.x  = element_text(size = 9),
          axis.text.y  = element_text(face = "plain"),
          panel.spacing.x = unit(0.2, "lines"))
  ggsave(outfile, p_main, height = 8, width = 12)
  cat("    wrote", basename(outfile), "\n")
}
plot_sen_decomp_grid(decomp_summary, cohort_n,
                     n_universe = median(ranked_meta$n_genes),
                     file.path(EXPLORATORY_PLOT_DIR, "senescence_decomposition_grid.pdf"))
plot_sen_decomp_rosette <- function(df, outfile) {
  cat_order_full <- c(
    "H1 Cell-cycle / arrest", "H2 DNA damage / DDR", "H3 Nuclear / chromatin",
    "H4 Anti-apoptotic / death receptors",
    "H5a SASP cytokines", "H5b SASP chemokines", "H5c SASP MMPs / TIMPs",
    "H5d SASP growth factors / IGFBPs", "H5e SASP serpins / fibrinolytic",
    "H5f SASP secreted misc",
    "H6 Lysosomal / SA-betagal", "H7 Metabolic / mitochondrial / ROS",
    "H8 Cell-surface markers", "ECM / matrix",
    "Non-can: Antiviral / interferon", "Non-can: Proteostasis / chaperone",
    "Non-can: WNT signaling", "Non-can: MAPK / stress signaling",
    "Non-can: Vascular / endothelial", "Non-can: Immune surface / adhesion",
    "Non-can: AD-specific", "Non-can: Skin / keratinization",
    "Non-can: Keratin tissue markers", "Non-can: Cadherin tissue markers",
    "Non-can: Pregnancy-specific", "Other (uncategorized)")
  pretty <- function(x) {
    x |> str_replace("\\u03B2", "beta") |>
      str_replace("^H[0-9][a-z]? ", "") |>
      str_replace("^Non-can: ", "") |>
      str_replace("^B[0-9]+ ", "")
  }
  label_lookup <- c(
    "H1 Cell-cycle / arrest" = "Cell cycle",
    "H2 DNA damage / DDR" = "DNA damage",
    "H3 Nuclear / chromatin" = "Chromatin",
    "H4 Anti-apoptotic / death receptors" = "Anti-apoptotic",
    "H5a SASP cytokines" = "SASP cytokines",
    "H5b SASP chemokines" = "SASP chemok.",
    "H5c SASP MMPs / TIMPs" = "SASP MMPs",
    "H5d SASP growth factors / IGFBPs" = "SASP GF/IGFBP",
    "H5e SASP serpins / fibrinolytic" = "SASP serpins",
    "H5f SASP secreted misc" = "SASP misc",
    "H6 Lysosomal / SA-betagal" = "Lysosomal",
    "H7 Metabolic / mitochondrial / ROS" = "Metabolic/ROS",
    "H8 Cell-surface markers" = "Cell surface",
    "ECM / matrix" = "ECM/matrix",
    "Non-can: Antiviral / interferon" = "Antiviral/IFN",
    "Non-can: Proteostasis / chaperone" = "Proteostasis",
    "Non-can: WNT signaling" = "WNT",
    "Non-can: MAPK / stress signaling" = "MAPK/stress",
    "Non-can: Vascular / endothelial" = "Vascular",
    "Non-can: Immune surface / adhesion" = "Immune adh.",
    "Non-can: AD-specific" = "AD-specific",
    "Non-can: Skin / keratinization" = "Skin/KRT",
    "Non-can: Keratin tissue markers" = "KRT markers",
    "Non-can: Cadherin tissue markers" = "Cadherin",
    "Non-can: Pregnancy-specific" = "Pregnancy",
    "Other (uncategorized)" = "Other"
  )
  key_df <- tibble(category_full = cat_order_full,
                   label = unname(label_lookup[cat_order_full]))

  d <- df |>
    mutate(category_full = as.character(category)) |>
    left_join(key_df, by = "category_full") |>
    mutate(category_pretty = factor(label, levels = key_df$label),
           Disease = comparison_disease(Comparison),
           Stratum = comparison_stratum(Comparison),
           Disease = factor(Disease, levels = DISEASE_ORDER),
           Stratum = factor(Stratum, levels = c("M", "F", "All")))
  cats <- key_df$label
  n_cats <- length(cats)
  d <- d |> mutate(x_num = as.integer(category_pretty),
                   y_num = as.integer(Stratum))
  labels_df <- key_df |>
    mutate(category_pretty = factor(label, levels = cats),
           x_num = seq_len(n_cats)) |>
    mutate(angle_raw = 90 - 360 * (x_num - 0.5) / n_cats,
           flip      = angle_raw < -90 | angle_raw > 90,
           angle     = ifelse(flip, angle_raw + 180, angle_raw),
           hjust     = ifelse(flip, 1, 0))

  outer_ring <- 3.5
  label_y    <- 4.62
  n_disease  <- n_distinct(na.omit(as.character(d$Disease)))
  facet_cols <- min(3, max(1, n_disease))

  p <- ggplot(d, aes(x = x_num, y = y_num, fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.3, width = 1, height = 1) +
    geom_text(data = labels_df,
              aes(x = x_num, y = label_y, label = label,
                  angle = angle, hjust = hjust),
              inherit.aes = FALSE, size = 3.1, color = "grey15") +
    facet_wrap(~ Disease, ncol = facet_cols) +
    coord_polar(theta = "x", clip = "off") +
    scale_x_continuous(limits = c(0.5, n_cats + 0.5), expand = c(0, 0),
                       breaks = NULL) +
    scale_y_continuous(limits = c(0, outer_ring + 1.95), expand = c(0, 0),
                       breaks = NULL) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-1, 1),
                         oob = scales::squish, name = "mean z") +
    labs(title = "Senescence decomposition by disease",
         caption = "Rings: outer = Combined | middle = Female | inner = Male",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title    = element_text(face = "bold"),
          strip.text    = element_text(face = "bold", size = 12),
          axis.text     = element_blank(),
          axis.ticks    = element_blank(),
          panel.grid    = element_blank(),
          legend.position = "right",
          plot.caption = element_text(size = 9, color = "grey25", hjust = 0),
          panel.spacing.x = unit(2.0, "lines"),
          panel.spacing.y = unit(1.5, "lines"),
          plot.margin = margin(8, 16, 4, 16))
  ggsave(outfile, p, height = ifelse(n_disease <= 3, 8.1, 11.6), width = 17.5)
  cat("    wrote", basename(outfile), "\n")
}
plot_sen_decomp_rosette(decomp_summary,
  file.path(PLOT_DIR, "senescence_decomposition_circular.pdf"))
plot_senescence_pathway_dotplot <- function(focused_df, outfile) {
  sen_pathways_themed <- tribble(
    ~pathway,                                                          ~theme,                       ~order,
    "Self-defined_Senescence",                                          "Custom + published refs",   1L,
    "SAUL_SEN_MAYO",                                                    "Custom + published refs",   2L,
    "REACTOME_CELLULAR_SENESCENCE",                                     "Senescence subtypes",       3L,
    "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",          "Senescence subtypes",       4L,
    "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE",                     "Senescence subtypes",       5L,
    "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE",           "Senescence subtypes",       6L,
    "REACTOME_ONCOGENE_INDUCED_SENESCENCE",                             "Senescence subtypes",       7L,
    "REACTOME_FORMATION_OF_SENESCENCE_ASSOCIATED_HETEROCHROMATIN_FOCI_SAHF","Senescence subtypes",   8L,
    "HALLMARK_DNA_REPAIR",                                              "H1/H2 cell-cycle / DDR",   10L,
    "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53",                       "H1/H2 cell-cycle / DDR",   11L,
    "REACTOME_REGULATION_OF_TP53_ACTIVITY",                              "H1/H2 cell-cycle / DDR",   12L,
    "REACTOME_TP53_REGULATES_METABOLIC_GENES",                           "H1/H2 cell-cycle / DDR",   13L,
    "HALLMARK_P53_PATHWAY",                                              "H1/H2 cell-cycle / DDR",   14L,
    "HALLMARK_APOPTOSIS",                                                "H4 anti-apoptotic",        20L,
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                                  "H5 SASP / inflammation",   30L,
    "HALLMARK_INFLAMMATORY_RESPONSE",                                    "H5 SASP / inflammation",   31L,
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",                                  "H5 SASP / inflammation",   32L,
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",                                "H7 metabolic / ROS",       40L,
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",                          "H7 metabolic / ROS",       41L,
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",                                "Aging-context",            50L,
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",                                "Aging-context",            51L,
    "HALLMARK_HYPOXIA",                                                  "Aging-context",            52L
  )
  theme_pretty <- function(t) {
    t |>
      str_replace("^H1/H2 ", "") |>
      str_replace("^H4 ", "") |>
      str_replace("^H5 ", "") |>
      str_replace("^H7 ", "") |>
      str_replace("^Aging-context", "Aging-context")
  }
  theme_levels_full <- c("Custom + published refs", "Senescence subtypes",
                          "H1/H2 cell-cycle / DDR", "H4 anti-apoptotic",
                          "H5 SASP / inflammation", "H7 metabolic / ROS",
                          "Aging-context")
  theme_levels_pretty <- theme_pretty(theme_levels_full)

  d <- focused_df |>
    inner_join(sen_pathways_themed, by = "pathway") |>
    mutate(
      theme_full = factor(theme, levels = theme_levels_full),
      theme = factor(theme_pretty(theme), levels = theme_levels_pretty),
      label = str_replace_all(pathway, "_", " ") |> str_trunc(60),
      label = factor(label,
                     levels = sen_pathways_themed |>
                       arrange(theme, order) |>
                       mutate(lbl = str_replace_all(pathway, "_", " ") |> str_trunc(60)) |>
                       pull(lbl)),
      Disease = comparison_disease(Comparison),
      Stratum = comparison_stratum(Comparison),
      Disease = factor(Disease, levels = DISEASE_ORDER),
      Stratum = factor(Stratum, levels = c("All","F","M")),
      sig = padj < 0.1)

  n_universe <- median(ranked_meta$n_genes)

  p_main <- ggplot(d, aes(x = Stratum, y = label,
                          fill = NES, size = -log10(pmax(padj, 1e-300)))) +
    geom_point(shape = 21, color = "grey20", stroke = 0.3) +
    geom_point(data = filter(d, sig), shape = 21, fill = NA,
               color = "black", stroke = 1.0) +
    facet_grid(theme ~ Disease, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-2.2, 2.2),
                         oob = scales::squish, name = "NES") +
    scale_size_continuous(range = c(2, 9), name = "-log10(padj)") +
    labs(title = "Senescence-focused pathway enrichment",
         x = NULL, y = NULL) +
    theme_clean +
    theme(strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
          strip.text.x      = element_text(face = "bold"),
          strip.placement   = "outside",
          axis.text.x       = element_text(size = 9),
          panel.spacing.y   = unit(0.4, "lines"),
          panel.spacing.x   = unit(0.3, "lines"))
  ggsave(outfile, p_main, height = 9, width = 12)
  cat("    wrote", basename(outfile), "\n")
}
plot_senescence_pathway_dotplot(focused_results,
  file.path(PLOT_DIR, "senescence_pathways_dotplot.pdf"))
cat("\n=== 7.6  Decomposing Self-defined_Bone by bone-biology sub-theme ===\n")

classify_bone_subtheme <- function(g) {
  if (g %in% c("RUNX2","SOX9","BGLAP","GCM2","IBSP","SPP1","MEPE","DMP1"))
    return("B1 Osteoblast TF / mineralization")
  if (g %in% c("BMP1","BMP2","BMP3","BMP4","BMP6","BMP7","BMPR1A",
               "ACVR1","NOG","GREM1","SMAD1","SMAD5"))
    return("B2 BMP signaling")
  if (g %in% c("WNT5A","WNT5B","WNT16","WNT7B","LRP5","LRP6","FZD9",
               "SOST","CTNNB1","SFRP1","WIF1","CCN4"))
    return("B3 WNT signaling")
  if (g %in% c("TGFB1","TGFBR1","SMAD2","SMAD3"))
    return("B4 TGF-beta signaling")
  if (g %in% c("NOTCH1","NOTCH2","NOTCH3","JAG1","DLL1"))
    return("B5 NOTCH signaling")
  if (g %in% c("CTSK","TNFSF11","TNFRSF11A","TNFRSF11B","ACP5","TRAP",
               "CALCR","CSF1"))
    return("B6 Bone resorption (osteoclast)")
  if (g %in% c("PTH","PTHLH","PTHD3","CALCR","CALCA","FGF23","ALPP","ALPL",
               "PHEX","ENPP1","CASR","DMP1"))
    return("B7 Endocrine / calcium-phosphate")
  if (str_detect(g, "^COL[0-9]") || g %in% c("PLOD2","CHST6","P4HB","LOX"))
    return("B8 Collagens & processing")
  if (str_detect(g, "^MMP[0-9]") || g %in% c("ADAM19","HPSE","TIMP3"))
    return("B9 MMPs / ECM remodeling")
  if (g %in% c("SPARC","POSTN","MGP","FBN1","FN1","THBS1","BSG","VTN",
               "HAS2","DCN","ANKRD1"))
    return("B10 Matrix / matricellular")
  if (str_detect(g, "^IGFBP[0-9]") ||
      g %in% c("IGF1","FGF2","FGF7","FGF19","VEGFA","EFNA2","PDGFA","NGF","HGF"))
    return("B11 IGFBPs / growth factors")
  if (g %in% c("IL1B","IL17A","IL6","TNF","S100A12","VCAM1","TNFSF8") ||
      str_detect(g, "^IL[0-9]"))
    return("B12 Inflammation / SASP overlap")
  if (g %in% c("ANXA1","ANXA2","ANXA5"))
    return("B13 Annexins / Ca-dependent membrane")
  if (g %in% c("CCN1","CCN2","CCN3","CCN4","CCN5","CCN6"))
    return("B14 CCN matricellular family")
  if (g %in% c("CD48","CD69","CD163","MS4A1","CEACAM8","CLEC12A","BSG","AHSP"))
    return("B15 Immune surface (bone-overlap)")
  if (g %in% c("CHI3L1","RGN","OLAH","OLR1","S100A12","NBR1","PTHD3",
               "JMJD6","UHRF1","KDM5D","ARID5B"))
    return("B16 Other bone-relevant secreted")
  return("Other (uncategorized)")
}
bone_custom <- gs_custom[["Self-defined_Bone"]]
bone_categories <- tibble(gene = bone_custom,
                          category = sapply(bone_custom, classify_bone_subtheme))

bone_decomp_rows <- list()
for (rname in names(ranked_lists)) {
  ranks <- ranked_lists[[rname]]
  d <- tibble(gene = names(ranks), z = as.numeric(ranks)) |>
    inner_join(bone_categories, by = "gene") |>
    mutate(Comparison = rname)
  bone_decomp_rows[[rname]] <- d
}
bone_decomp_df <- bind_rows(bone_decomp_rows)

bone_decomp_summary <- bone_decomp_df |>
  group_by(Comparison, category) |>
  summarise(n_in_assay = n(),
            mean_z     = mean(z),
            median_z   = median(z),
            n_up       = sum(z > 0),
            n_down     = sum(z < 0),
            pct_up     = 100 * sum(z > 0) / n(),
            .groups = "drop") |>
  arrange(category, Comparison)

write_csv(bone_decomp_summary,
          file.path(OUT_DIR, "bone_decomposition.csv"))
cat(sprintf("  wrote bone_decomposition.csv (%d rows)\n",
            nrow(bone_decomp_summary)))
plot_bone_decomp <- function(df, outfile) {
  cat_order <- c("B1 Osteoblast TF / mineralization",
                 "B2 BMP signaling", "B3 WNT signaling", "B4 TGF-beta signaling",
                 "B5 NOTCH signaling",
                 "B6 Bone resorption (osteoclast)",
                 "B7 Endocrine / calcium-phosphate",
                 "B8 Collagens & processing",
                 "B9 MMPs / ECM remodeling",
                 "B10 Matrix / matricellular",
                 "B11 IGFBPs / growth factors",
                 "B12 Inflammation / SASP overlap",
                 "B13 Annexins / Ca-dependent membrane",
                 "B14 CCN matricellular family",
                 "B15 Immune surface (bone-overlap)",
                 "B16 Other bone-relevant secreted",
                 "Other (uncategorized)")
  d <- df |> mutate(category = factor(category, levels = rev(cat_order)),
                    label = sprintf("%s (n=%d)", category, n_in_assay))
  if (nrow(d) == 0) return(invisible())
  p <- ggplot(d, aes(x = Comparison, y = category, fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.2f\n(n=%d)", mean_z, n_in_assay)),
              size = 3, color = "black") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-0.6, 0.6),
                         oob = scales::squish, name = "mean z") +
    labs(title = "Decomposition of Self-defined bone panel by bone-biology sub-theme",
         subtitle = "Mean z-statistic across genes in each sub-theme. Cells show mean z and n genes measured.",
         x = NULL, y = NULL) +
    theme_clean +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(face = "bold"))
  ggsave(outfile, p, height = 7, width = 9)
  cat("    wrote", basename(outfile), "\n")
}
plot_bone_decomp(bone_decomp_summary,
  file.path(EXPLORATORY_PLOT_DIR, "bone_decomposition.pdf"))
plot_bone_decomp_grid <- function(df, cohort_n, n_universe, outfile) {
  cat_order_full <- c("B1 Osteoblast TF / mineralization",
                      "B2 BMP signaling", "B3 WNT signaling", "B4 TGF-beta signaling",
                      "B5 NOTCH signaling",
                      "B6 Bone resorption (osteoclast)",
                      "B7 Endocrine / calcium-phosphate",
                      "B8 Collagens & processing", "B9 MMPs / ECM remodeling",
                      "B10 Matrix / matricellular",
                      "B11 IGFBPs / growth factors",
                      "B12 Inflammation / SASP overlap",
                      "B13 Annexins / Ca-dependent membrane",
                      "B14 CCN matricellular family",
                      "B15 Immune surface (bone-overlap)",
                      "B16 Other bone-relevant secreted",
                      "Other (uncategorized)")
  pretty <- function(x) x |> str_replace("\\u03B2", "beta") |> str_replace("^B[0-9]+ ", "")
  d <- df |>
    mutate(category_full = str_replace(category, "\\u03B2", "beta"),
           category_full = factor(category_full, levels = cat_order_full),
           category_pretty = factor(pretty(as.character(category_full)),
                                    levels = unique(pretty(cat_order_full))),
           Disease = comparison_disease(Comparison),
           Stratum = comparison_stratum(Comparison),
           Disease = factor(Disease, levels = DISEASE_ORDER),
           Stratum = factor(Stratum, levels = c("All","F","M")))

  p_main <- ggplot(d, aes(x = Stratum, y = forcats::fct_rev(category_pretty),
                          fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%+.2f", mean_z)),
              color = "grey15", size = 2.4) +
    facet_grid(. ~ Disease, switch = "y") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-1, 1),
                         oob = scales::squish, name = "mean z") +
    labs(title = "Self-defined bone panel - gene-level decomposition by bone-biology sub-theme",
         x = NULL, y = NULL) +
    theme_clean +
    theme(strip.text.x = element_text(face = "bold"),
          axis.text.x  = element_text(size = 9),
          axis.text.y  = element_text(face = "plain"),
          panel.spacing.x = unit(0.2, "lines"))
  ggsave(outfile, p_main, height = 7, width = 12)
  cat("    wrote", basename(outfile), "\n")
}
plot_bone_decomp_grid(bone_decomp_summary, cohort_n,
                      n_universe = median(ranked_meta$n_genes),
                      file.path(EXPLORATORY_PLOT_DIR, "bone_decomposition_grid.pdf"))
plot_bone_decomp_rosette <- function(df, outfile) {
  cat_order_full <- c("B1 Osteoblast TF / mineralization",
                      "B2 BMP signaling", "B3 WNT signaling", "B4 TGF-beta signaling",
                      "B5 NOTCH signaling",
                      "B6 Bone resorption (osteoclast)",
                      "B7 Endocrine / calcium-phosphate",
                      "B8 Collagens & processing", "B9 MMPs / ECM remodeling",
                      "B10 Matrix / matricellular",
                      "B11 IGFBPs / growth factors",
                      "B12 Inflammation / SASP overlap",
                      "B13 Annexins / Ca-dependent membrane",
                      "B14 CCN matricellular family",
                      "B15 Immune surface (bone-overlap)",
                      "B16 Other bone-relevant secreted",
                      "Other (uncategorized)")
  pretty <- function(x) x |> str_replace("\\u03B2", "beta") |> str_replace("^B[0-9]+ ", "")
  label_lookup <- c(
    "B1 Osteoblast TF / mineralization" = "Osteoblast TF",
    "B2 BMP signaling" = "BMP",
    "B3 WNT signaling" = "WNT",
    "B4 TGF-beta signaling" = "TGF-beta",
    "B5 NOTCH signaling" = "NOTCH",
    "B6 Bone resorption (osteoclast)" = "Resorption",
    "B7 Endocrine / calcium-phosphate" = "Ca/phosphate",
    "B8 Collagens & processing" = "Collagens",
    "B9 MMPs / ECM remodeling" = "MMP/ECM",
    "B10 Matrix / matricellular" = "Matrix",
    "B11 IGFBPs / growth factors" = "IGFBP/GF",
    "B12 Inflammation / SASP overlap" = "Inflamm/SASP",
    "B13 Annexins / Ca-dependent membrane" = "Annexins/Ca",
    "B14 CCN matricellular family" = "CCN",
    "B15 Immune surface (bone-overlap)" = "Immune surface",
    "B16 Other bone-relevant secreted" = "Other secreted",
    "Other (uncategorized)" = "Other"
  )
  key_df <- tibble(category_full = cat_order_full,
                   label = unname(label_lookup[cat_order_full]))

  d <- df |>
    mutate(category_full = str_replace(category, "\\u03B2", "beta"),
           category_full = as.character(category_full)) |>
    left_join(key_df, by = "category_full") |>
    mutate(category_pretty = factor(label, levels = key_df$label),
           Disease = comparison_disease(Comparison),
           Stratum = comparison_stratum(Comparison),
           Disease = factor(Disease, levels = DISEASE_ORDER),
           Stratum = factor(Stratum, levels = c("M", "F", "All")))

  cats <- key_df$label
  n_cats <- length(cats)
  d <- d |> mutate(x_num = as.integer(category_pretty),
                   y_num = as.integer(Stratum))
  labels_df <- key_df |>
    mutate(category_pretty = factor(label, levels = cats),
           x_num = seq_len(n_cats)) |>
    mutate(angle_raw = 90 - 360 * (x_num - 0.5) / n_cats,
           flip      = angle_raw < -90 | angle_raw > 90,
           angle     = ifelse(flip, angle_raw + 180, angle_raw),
           hjust     = ifelse(flip, 1, 0))

  outer_ring <- 3.5
  label_y    <- 4.58
  n_disease  <- n_distinct(na.omit(as.character(d$Disease)))
  facet_cols <- min(3, max(1, n_disease))

  p <- ggplot(d, aes(x = x_num, y = y_num, fill = mean_z)) +
    geom_tile(color = "white", linewidth = 0.3, width = 1, height = 1) +
    geom_text(data = labels_df,
              aes(x = x_num, y = label_y, label = label,
                  angle = angle, hjust = hjust),
              inherit.aes = FALSE, size = 3.35, color = "grey15") +
    facet_wrap(~ Disease, ncol = facet_cols) +
    coord_polar(theta = "x", clip = "off") +
    scale_x_continuous(limits = c(0.5, n_cats + 0.5), expand = c(0, 0),
                       breaks = NULL) +
    scale_y_continuous(limits = c(0, outer_ring + 1.85), expand = c(0, 0),
                       breaks = NULL) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-1, 1),
                         oob = scales::squish, name = "mean z") +
    labs(title = "Bone-biology decomposition by disease",
         caption = "Rings: outer = Combined | middle = Female | inner = Male",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title    = element_text(face = "bold"),
          strip.text    = element_text(face = "bold", size = 12),
          axis.text     = element_blank(),
          axis.ticks    = element_blank(),
          panel.grid    = element_blank(),
          legend.position = "right",
          plot.caption = element_text(size = 9, color = "grey25", hjust = 0),
          panel.spacing.x = unit(2.0, "lines"),
          panel.spacing.y = unit(1.5, "lines"),
          plot.margin = margin(8, 16, 4, 16))
  ggsave(outfile, p, height = ifelse(n_disease <= 3, 7.7, 11.0), width = 17.5)
  cat("    wrote", basename(outfile), "\n")
}
plot_bone_decomp_rosette(bone_decomp_summary,
  file.path(PLOT_DIR, "bone_decomposition_circular.pdf"))
plot_bone_pathway_dotplot <- function(focused_df, outfile) {
  bone_pathways_themed <- tribble(
    ~pathway,                                                 ~theme,                          ~order,
    "Self-defined_Bone",                                      "Custom panel",                  1L,
    "REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT",              "B1 RUNX2 / osteoblast",         2L,
    "REACTOME_RUNX2_REGULATES_OSTEOBLAST_DIFFERENTIATION",    "B1 RUNX2 / osteoblast",         3L,
    "REACTOME_SIGNALING_BY_BMP",                              "B2 BMP signaling",              4L,
    "REACTOME_SIGNALING_BY_WNT",                              "B3 WNT signaling",              5L,
    "REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT",    "B3 WNT signaling",              6L,
    "REACTOME_SIGNALING_BY_WNT_IN_CANCER",                    "B3 WNT signaling",              7L,
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",                    "B3 WNT signaling",              8L,
    "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",        "B4 TGF-beta signaling",            9L,
    "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",   "B4 TGF-beta signaling",           10L,
    "HALLMARK_TGF_BETA_SIGNALING",                            "B4 TGF-beta signaling",           11L,
    "REACTOME_SIGNALING_BY_NOTCH",                            "B5 NOTCH signaling",           12L,
    "HALLMARK_NOTCH_SIGNALING",                               "B5 NOTCH signaling",           13L,
    "REACTOME_COLLAGEN_FORMATION",                            "B8 Collagens",                 14L,
    "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",   "B8 Collagens",                 15L,
    "REACTOME_COLLAGEN_DEGRADATION",                          "B8 Collagens",                 16L,
    "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",                  "B8 Collagens",                 17L,
    "REACTOME_ECM_PROTEOGLYCANS",                             "B10 ECM / integrin",           18L,
    "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",            "B10 ECM / integrin",           19L,
    "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",        "B10 ECM / integrin",           20L,
    "REACTOME_CHONDROITIN_SULFATE_BIOSYNTHESIS",              "B10 ECM / integrin",           21L,
    "REACTOME_CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM","B10 ECM / integrin",          22L,
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",                       "Sex hormones",                 30L,
    "HALLMARK_ESTROGEN_RESPONSE_LATE",                        "Sex hormones",                 31L,
    "HALLMARK_ANDROGEN_RESPONSE",                             "Sex hormones",                 32L
  )
  theme_pretty <- function(t) t |> str_replace("^B[0-9]+ ", "") |>
    str_replace("\\u03B2", "beta")
  theme_levels_full <- c("Custom panel",
                          "B1 RUNX2 / osteoblast", "B2 BMP signaling",
                          "B3 WNT signaling", "B4 TGF-beta signaling",
                          "B5 NOTCH signaling", "B8 Collagens",
                          "B10 ECM / integrin", "Sex hormones")
  theme_levels_pretty <- theme_pretty(theme_levels_full)

  d <- focused_df |>
    inner_join(bone_pathways_themed, by = "pathway") |>
    mutate(
      theme_full = theme,
      theme_full = ifelse(theme_full == "B4 TGF-beta signaling", "B4 TGF-beta signaling", theme_full),
      theme = theme_pretty(theme_full),
      theme = factor(theme, levels = theme_levels_pretty),
      label = str_replace_all(pathway, "_", " ") |> str_trunc(60),
      label = factor(label,
                     levels = bone_pathways_themed |>
                       arrange(theme, order) |>
                       mutate(lbl = str_replace_all(pathway, "_", " ") |> str_trunc(60)) |>
                       pull(lbl)),
      Disease = comparison_disease(Comparison),
      Stratum = comparison_stratum(Comparison),
      Disease = factor(Disease, levels = DISEASE_ORDER),
      Stratum = factor(Stratum, levels = c("All","F","M")),
      sig = padj < 0.1)

  n_universe <- median(ranked_meta$n_genes)

  p_main <- ggplot(d, aes(x = Stratum, y = label,
                          fill = NES, size = -log10(pmax(padj, 1e-300)))) +
    geom_point(shape = 21, color = "grey20", stroke = 0.3) +
    geom_point(data = filter(d, sig), shape = 21, fill = NA,
               color = "black", stroke = 1.0) +
    facet_grid(theme ~ Disease, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1432b",
                         midpoint = 0, limits = c(-2.2, 2.2),
                         oob = scales::squish, name = "NES") +
    scale_size_continuous(range = c(2, 9), name = "-log10(padj)") +
    labs(title = "Bone-focused pathway enrichment",
         x = NULL, y = NULL) +
    theme_clean +
    theme(strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
          strip.text.x      = element_text(face = "bold"),
          strip.placement   = "outside",
          axis.text.x       = element_text(size = 9),
          panel.spacing.y   = unit(0.4, "lines"),
          panel.spacing.x   = unit(0.3, "lines"))
  ggsave(outfile, p_main, height = 9, width = 12)
  cat("    wrote", basename(outfile), "\n")
}
plot_bone_pathway_dotplot(focused_results,
  file.path(PLOT_DIR, "bone_pathways_dotplot.pdf"))
cat("\n=== 8. Summary ===\n")

cat("\n--- Tier 1 (PRIMARY) significant pathways per stratum ---\n")
focused_sig |>
  count(Comparison, Source, name = "n_sig") |>
  pivot_wider(names_from = Source, values_from = n_sig, values_fill = 0) |>
  print(n = Inf)

if (nrow(focused_sig) > 0) {
  cat("\n--- Tier 1 hits (all, sorted by padj) ---\n")
  focused_sig |>
    select(Comparison, Source, Theme = Source, pathway, NES, padj, size) |>
    print(n = Inf)
}

cat("\n--- Custom panel rank in Tier 1 vs Tier 2 ---\n")
position_table |>
  select(Universe, Comparison, pathway, NES, pval, padj, rank, total, percentile) |>
  print(n = Inf)

cat("\n--- Tier 2 (SUPPLEMENTARY) sig counts ---\n")
unified_sig |>
  count(Comparison, Source, name = "n_sig") |>
  pivot_wider(names_from = Source, values_from = n_sig, values_fill = 0) |>
  print(n = Inf)

cat("\nDone. Outputs are under:\n  ", OUT_DIR, "\n")
