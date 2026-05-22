#!/usr/bin/env Rscript

Sys.setenv(VROOM_CONNECTION_SIZE = 2^24)

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(variancePartition)
  library(BiocParallel)
  library(limma)
  library(yaml)
})

script_dir <- {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    normalizePath(getwd())
  }
}

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

as_null <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  if (is.na(x) || identical(x, "")) return(NULL)
  if (toupper(as.character(x)) %in% c("NA", "NULL", "NONE")) return(NULL)
  x
}

as_logical <- function(x, default = FALSE) {
  x <- as_null(x)
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x)
}

as_numeric_or_null <- function(x) {
  x <- as_null(x)
  if (is.null(x)) return(NULL)
  as.numeric(x)
}

required_config <- function(x, label) {
  x <- as_null(x)
  if (is.null(x)) stop("Set ", label, " in the analysis config.")
  as.character(x)
}

read_private_ids <- function(path) {
  if (is.null(path) || !file.exists(path)) return(character())
  ids <- trimws(readLines(path, warn = FALSE))
  ids[nzchar(ids)]
}

deep_modify <- function(x, y) {
  if (length(y) == 0L) return(x)
  for (nm in names(y)) {
    if (is.list(x[[nm]]) && is.list(y[[nm]])) {
      x[[nm]] <- deep_modify(x[[nm]], y[[nm]])
    } else {
      x[nm] <- list(y[[nm]])
    }
  }
  x
}

default_config <- list(
  paths = list(
    root_dir = ".",
    data_dir = "data",
    output_base = "outputs",
    expression_file = "sample_protein_expression.csv",
    sample_metadata_file = "sample_metadata.csv",
    demographics_file = "sample_demographics.xlsx",
    protein_metadata_file = "protein_metadata.csv",
    panel_file = Sys.getenv("MANUSCRIPT_PANEL_SOMAKEYS", unset = ""),
    output_dir = NULL
  ),
  analysis = list(
    group_mode = "individual",
    min_baseline_age = 55,
    sex_filter = "all",
    subject_id_field = NULL,
    exclude_subject_ids_file = "private/exclude_subject_ids.txt",
    exclude_subject_ids = NULL,
    save_fit = TRUE,
    n_cores = 1
  ),
  plots = list(
    run_full_volcano = TRUE,
    run_panel_volcano = TRUE,
    p_threshold = 0.05,
    log2fc_guide = log2(1.1),
    color_higher = "#D55E00",
    color_lower = "#0072B2",
    color_neutral = "#a6adb7"
  ),
  test = list(
    test_n_somakeys = NULL,
    test_panel_somakeys_per_panel = 20
  )
)

args <- parse_args()
config_file <- if (!is.null(args$config)) {
  args$config
} else {
  candidate <- file.path(script_dir, "config.yml")
  if (file.exists(candidate)) candidate else file.path(script_dir, "config_template.yml")
}

cfg <- default_config
if (file.exists(config_file)) {
  cfg <- deep_modify(cfg, yaml::read_yaml(config_file))
} else {
  warning("Config file not found; using built-in defaults: ", config_file)
}

for (nm in setdiff(names(args), "config")) {
  val <- args[[nm]]
  if (nm %in% names(cfg$paths)) cfg$paths[nm] <- list(val)
  if (nm %in% names(cfg$analysis)) cfg$analysis[nm] <- list(val)
  if (nm %in% names(cfg$plots)) cfg$plots[nm] <- list(val)
  if (nm %in% names(cfg$test)) cfg$test[nm] <- list(val)
}

root_dir <- normalizePath(cfg$paths$root_dir, mustWork = FALSE)
resolve_path <- function(base, path) {
  path <- as_null(path)
  if (is.null(path)) return(NULL)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) normalizePath(path, mustWork = FALSE)
  else normalizePath(file.path(base, path), mustWork = FALSE)
}

data_dir <- resolve_path(root_dir, cfg$paths$data_dir)
output_base <- resolve_path(root_dir, cfg$paths$output_base)
expr_file <- file.path(data_dir, cfg$paths$expression_file)
smeta_file <- file.path(data_dir, cfg$paths$sample_metadata_file)
demo_file <- file.path(data_dir, cfg$paths$demographics_file)
protein_file <- file.path(data_dir, cfg$paths$protein_metadata_file)
panel_file_value <- as_null(cfg$paths$panel_file)
if (identical(panel_file_value, "/path/to/manuscript_panel_somakeys.csv")) {
  panel_file_value <- NULL
}
if (is.null(panel_file_value)) {
  panel_file_value <- as_null(Sys.getenv("MANUSCRIPT_PANEL_SOMAKEYS", unset = ""))
}
panel_file <- resolve_path(root_dir, panel_file_value)

group_mode <- as.character(cfg$analysis$group_mode)
if (!group_mode %in% c("combined", "individual")) {
  stop("analysis.group_mode must be 'combined' or 'individual'.")
}

min_baseline_age <- as_numeric_or_null(cfg$analysis$min_baseline_age)
sex_filter <- switch(
  tolower(as.character(cfg$analysis$sex_filter)),
  "all" = "all",
  "none" = "all",
  "f" = "F",
  "female" = "F",
  "m" = "M",
  "male" = "M",
  stop("analysis.sex_filter must be 'all', 'F', or 'M'.")
)
sex_filter_label <- switch(sex_filter, "all" = "all", "F" = "female", "M" = "male")
subject_id_field <- required_config(cfg$analysis$subject_id_field, "analysis.subject_id_field")
exclude_subject_ids_file <- resolve_path(root_dir, cfg$analysis$exclude_subject_ids_file)
exclude_subject_ids <- unique(c(
  as.character(unlist(cfg$analysis$exclude_subject_ids)),
  read_private_ids(exclude_subject_ids_file)
))
exclude_subject_ids <- exclude_subject_ids[!is.na(exclude_subject_ids) & nzchar(exclude_subject_ids)]
save_fit <- as_logical(cfg$analysis$save_fit, TRUE)
n_cores <- as.integer(cfg$analysis$n_cores)
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

run_full_volcano <- as_logical(cfg$plots$run_full_volcano, TRUE)
run_panel_volcano <- as_logical(cfg$plots$run_panel_volcano, TRUE)
p_threshold <- as.numeric(cfg$plots$p_threshold)
log2fc_guide <- as.numeric(cfg$plots$log2fc_guide)
test_n_somakeys <- as_numeric_or_null(cfg$test$test_n_somakeys)
test_panel_somakeys_per_panel <- as.integer(cfg$test$test_panel_somakeys_per_panel)
if (is.na(test_panel_somakeys_per_panel)) test_panel_somakeys_per_panel <- 0L

age_suffix <- if (is.null(min_baseline_age)) "" else sprintf("_baseAge%g", min_baseline_age)
sex_suffix <- if (sex_filter == "all") "" else sprintf("_sex%s", sex_filter)
default_out_dir <- file.path(output_base, sprintf("dream_%s%s%s", group_mode, age_suffix, sex_suffix))
out_dir <- resolve_path(root_dir, cfg$paths$output_dir)
if (is.null(out_dir)) out_dir <- default_out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output dir: ", out_dir)
message("Model group mode: ", group_mode)
message("Sex filter: ", sex_filter_label)

if (n_cores == 1L) {
  bp_param <- SerialParam()
  message("Using SerialParam for dream fitting.")
} else {
  bp_param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
  message("Using MulticoreParam with ", n_cores, " workers for dream fitting.")
}

message("[1/6] Loading data")
expr_raw <- read_csv(expr_file, show_col_types = FALSE)
smeta <- read_csv(smeta_file, show_col_types = FALSE, guess_max = 5000, name_repair = "unique")
demo <- read_excel(demo_file)
protein_dict <- read_csv(protein_file, show_col_types = FALSE)

sid_cols <- grep("^SampleId", names(smeta), value = TRUE)
if (length(sid_cols) > 1L) {
  stopifnot(identical(as.character(smeta[[sid_cols[1]]]), as.character(smeta[[sid_cols[2]]])))
  smeta[[sid_cols[2]]] <- NULL
  names(smeta)[names(smeta) == sid_cols[1]] <- "SampleId"
}

expr_raw$SampleId <- as.character(expr_raw$SampleId)
smeta$SampleId <- as.character(smeta$SampleId)
demo$SampleId <- as.character(demo$SampleId)

message("[2/6] QC filtering and log2 replicate averaging")
soma_cols <- grep("^SL[0-9]", names(expr_raw), value = TRUE)
hce_cols <- grep("^HCE", names(expr_raw), value = TRUE)
message("Detected ", length(soma_cols), " SOMAmer columns; ignored ", length(hce_cols), " HCE controls.")

if (!is.null(test_n_somakeys)) {
  test_n_somakeys <- min(as.integer(test_n_somakeys), length(soma_cols))
  test_keys <- head(soma_cols, test_n_somakeys)
  if (!is.null(panel_file) && file.exists(panel_file) &&
      test_panel_somakeys_per_panel > 0L) {
    panel_keys <- read_csv(panel_file, show_col_types = FALSE) %>%
      group_by(panel) %>%
      slice_head(n = test_panel_somakeys_per_panel) %>%
      ungroup() %>%
      pull(SomaKey)
    test_keys <- unique(c(test_keys, intersect(panel_keys, soma_cols)))
  }
  soma_cols <- test_keys
  message("TEST MODE: fitting ", length(soma_cols), " SOMAmer columns.")
}

keep_ids <- intersect(
  smeta$SampleId[smeta$RowCheck == "PASS"],
  smeta$SampleId[smeta$SampleType == "Sample"]
)
expr <- expr_raw %>% filter(SampleId %in% keep_ids)

soma_mat_raw <- as.matrix(expr[, soma_cols, drop = FALSE])
rownames(soma_mat_raw) <- expr$SampleId
stopifnot(all(!is.na(soma_mat_raw)), all(soma_mat_raw > 0))

log2_long <- log2(soma_mat_raw)
log2_byid <- rowsum(log2_long, group = rownames(log2_long), reorder = TRUE)
counts_byid <- as.numeric(table(rownames(log2_long))[rownames(log2_byid)])
expr_mat <- t(log2_byid / counts_byid)
message("Expression matrix after averaging: ", nrow(expr_mat), " proteins x ", ncol(expr_mat), " samples")

message("[3/6] Building covariate table")
if (!subject_id_field %in% names(demo)) {
  stop("Missing subject ID column in demographics file: ", subject_id_field)
}
demo_u <- demo %>%
  filter(SampleId %in% colnames(expr_mat)) %>%
  distinct(SampleId, .keep_all = TRUE) %>%
  mutate(
    source_subject_id = as.character(.data[[subject_id_field]]),
    Visit_int = as.integer(Visit)
  )

smeta_u <- smeta %>%
  filter(SampleId %in% colnames(expr_mat)) %>%
  distinct(SampleId, .keep_all = TRUE) %>%
  select(SampleId, SampleType, PlateId, PlateRunDate, ScannerID, RowCheck)

covar <- demo_u %>%
  left_join(smeta_u, by = "SampleId") %>%
  filter(!(source_subject_id %in% exclude_subject_ids)) %>%
  mutate(
    group_combined = case_when(
      Diagnosis_group == "HC" ~ "HC",
      Diagnosis_group %in% c("MCI", "AD") ~ "AD_spectrum",
      Diagnosis_group == "PD" ~ "PD",
      Diagnosis_group %in% c("MCI-PD", "LBD") ~ "LB_spectrum",
      TRUE ~ NA_character_
    ),
    group_combined = factor(group_combined, levels = c("HC", "AD_spectrum", "PD", "LB_spectrum")),
    group_individual = case_when(
      Diagnosis_group %in% c("HC", "MCI", "AD", "PD", "LBD") ~ Diagnosis_group,
      Diagnosis_group == "MCI-PD" ~ "MCI_PD",
      TRUE ~ NA_character_
    ),
    group_individual = factor(group_individual, levels = c("HC", "MCI", "AD", "PD", "MCI_PD", "LBD")),
    Sex = factor(Sex),
    SubjectID = factor(paste0("S", match(source_subject_id, unique(source_subject_id)))),
    age_at_visit = as.numeric(Age_years)
  ) %>%
  select(-source_subject_id)

if (length(exclude_subject_ids) > 0L) {
  message("Applied prespecified subject exclusions: ", length(exclude_subject_ids), " subject(s)")
}

if (sex_filter != "all") {
  covar <- covar %>% filter(Sex == sex_filter)
}

group_col <- switch(group_mode, "combined" = "group_combined", "individual" = "group_individual")
covar$group_var <- covar[[group_col]]

covar <- covar %>%
  filter(!is.na(group_var), !is.na(age_at_visit), !is.na(Sex))

if (!is.null(min_baseline_age)) {
  baseline_age <- covar %>%
    group_by(SubjectID) %>%
    summarise(baseline_age = min(age_at_visit, na.rm = TRUE), .groups = "drop")
  keep_subj <- baseline_age$SubjectID[baseline_age$baseline_age >= min_baseline_age]
  covar <- covar %>% filter(SubjectID %in% keep_subj)
}

covar$group_var <- droplevels(covar$group_var)
covar$Sex <- droplevels(covar$Sex)
covar$SubjectID <- droplevels(covar$SubjectID)

if (!"HC" %in% levels(covar$group_var)) {
  stop("HC group has no samples after filtering.")
}
if (length(setdiff(levels(covar$group_var), "HC")) == 0L) {
  stop("No disease group remains after filtering.")
}

expr_mat <- expr_mat[, intersect(colnames(expr_mat), covar$SampleId), drop = FALSE]
covar <- covar[match(colnames(expr_mat), covar$SampleId), ]
covar <- as.data.frame(covar)
rownames(covar) <- covar$SampleId
stopifnot(identical(colnames(expr_mat), covar$SampleId))

message("Final analysis set: ", nrow(expr_mat), " proteins x ", ncol(expr_mat),
        " samples / ", length(unique(covar$SubjectID)), " subjects")

group_sizes <- covar %>%
  group_by(group_var) %>%
  summarise(
    n_samples = n(),
    n_participants = n_distinct(SubjectID),
    .groups = "drop"
  ) %>%
  arrange(match(group_var, levels(covar$group_var)))
write_csv(group_sizes, file.path(out_dir, "group_sizes.csv"))

message("[4/6] Fitting dream LMM")
fixed_terms <- c("group_var", "age_at_visit")
if (sex_filter == "all" && nlevels(covar$Sex) >= 2L) {
  fixed_terms <- c(fixed_terms, "Sex")
}

form_txt <- paste("~", paste(fixed_terms, collapse = " + "))
if (any(table(covar$SubjectID) > 1L)) {
  form_txt <- paste(form_txt, "+ (1 | SubjectID)")
}
form <- as.formula(form_txt)
message("Model formula: ", deparse(form))

fit <- dream(exprObj = expr_mat, formula = form, data = covar, BPPARAM = bp_param)
fit <- eBayes(fit)

if (save_fit) {
  saveRDS(fit, file.path(out_dir, "dream_fit.rds"))
}
saveRDS(covar, file.path(out_dir, "covar_used.rds"))

message("[5/6] Writing disease-vs-HC top tables")
disease_levels <- setdiff(levels(covar$group_var), "HC")
coefs <- setNames(paste0("group_var", disease_levels), paste0(disease_levels, "_vs_HC"))
missing <- setdiff(coefs, colnames(coef(fit)))
if (length(missing) > 0L) {
  stop("Missing model coefficients: ", paste(missing, collapse = ", "))
}

ann <- protein_dict %>%
  select(SomaKey, SomaId, Target, TargetFullName, UniProt, EntrezGeneSymbol)

size_lookup <- setNames(group_sizes$n_samples, as.character(group_sizes$group_var))
model_label <- sprintf(
  "dream LMM; group_mode=%s; baseline_age>=%s; sex_stratum=%s; outcome=log2(RFU)",
  group_mode,
  if (is.null(min_baseline_age)) "none" else sprintf("%g", min_baseline_age),
  sex_filter_label
)

for (nm in names(coefs)) {
  tt <- topTable(fit, coef = coefs[[nm]], number = Inf, adjust.method = "BH", sort.by = "P")
  tt_ann <- tt %>%
    rownames_to_column("SomaKey") %>%
    left_join(ann, by = "SomaKey")

  disease_level <- sub("_vs_HC$", "", nm)
  tt_ann$contrast <- nm
  tt_ann$comparison_group <- disease_level
  tt_ann$reference_group <- "HC"
  tt_ann$group_mode <- group_mode
  tt_ann$min_baseline_age <- if (is.null(min_baseline_age)) NA_real_ else min_baseline_age
  tt_ann$sex_stratum <- sex_filter_label
  tt_ann$model <- model_label
  tt_ann$model_formula <- deparse(form)
  tt_ann$coef_name <- coefs[[nm]]
  tt_ann[[paste0("n_", disease_level)]] <- size_lookup[[disease_level]]
  tt_ann$n_HC <- size_lookup[["HC"]]

  write_csv(tt_ann, file.path(out_dir, sprintf("toptable_%s.csv", nm)))
  message("  ", nm, ": ", sum(tt_ann$adj.P.Val < 0.05, na.rm = TRUE), " proteins at FDR<0.05")
}

plate_bal <- table(covar$group_var, covar$PlateId, useNA = "ifany")
write.csv(addmargins(plate_bal), file.path(out_dir, "plate_x_group_balance.csv"))

message("[6/6] Plotting")
if (run_full_volcano) {
  plot_env <- new.env(parent = globalenv())
  plot_env$DREAM_OUT_DIR <- out_dir
  plot_env$PLOT_P_THRESHOLD <- p_threshold
  plot_env$PLOT_LOG2FC_GUIDE <- log2fc_guide
  plot_env$PLOT_COLOR_HIGHER <- cfg$plots$color_higher
  plot_env$PLOT_COLOR_LOWER <- cfg$plots$color_lower
  plot_env$PLOT_COLOR_NEUTRAL <- cfg$plots$color_neutral
  source(file.path(script_dir, "02_make_proteome_volcano.R"), local = plot_env)
}

if (run_panel_volcano && group_mode == "individual" &&
    !is.null(panel_file) && file.exists(panel_file)) {
  panel_env <- new.env(parent = globalenv())
  panel_env$DREAM_OUT_DIR <- out_dir
  panel_env$PANEL_SOURCE <- panel_file
  panel_env$PLOT_P_THRESHOLD <- p_threshold
  panel_env$PLOT_LOG2FC_GUIDE <- log2fc_guide
  panel_env$PLOT_COLOR_HIGHER <- cfg$plots$color_higher
  panel_env$PLOT_COLOR_LOWER <- cfg$plots$color_lower
  panel_env$PLOT_COLOR_NEUTRAL <- cfg$plots$color_neutral
  source(file.path(script_dir, "03_make_panel_volcano.R"), local = panel_env)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Results written to: ", out_dir)
