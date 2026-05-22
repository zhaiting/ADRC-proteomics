# Final Figure Code

Analysis code for the manuscript figures: plasma SomaScan proteomic association with neurodegenerative disease groups, correlation of bone- and senescence-panel proteins with plasma Alzheimer's disease biomarkers, and preranked pathway enrichment.

## Folders

- `proteome_association/`: fits protein association models and makes volcano plots from SomaScan data.
- `correlation/`: builds protein-biomarker correlation models and final heatmaps from the proteome association outputs.
- `pathway_analysis/`: runs preranked GSEA and pathway summaries from association statistics, using the bone and senescence gene lists provided in the manuscript.

## Inputs To Provide Locally

Create `data/` with the SomaScan expression, sample metadata, demographics, and protein annotation files. Create `biomarker_data/` with the plasma biomarker spreadsheets used by `correlation/correlation_model_setup.Rmd`. Copy `proteome_association/config_template.yml` to ignored `proteome_association/config.yml` and edit local file names and `analysis.subject_id_field`.

Provide the manuscript panel SomaKey table with `MANUSCRIPT_PANEL_SOMAKEYS`. The proteome association script uses this value unless `paths.panel_file` is set in `config.yml`.

Before running the correlation setup, set local column and file names:

``` bash
export SUBJECT_ID_FIELD="local_subject_id_column"
export BIOMARKER_ID_FIELD="local_biomarker_id_column"
export CORRELATION_EXPRESSION_FILE="sample_protein_expression.csv"
export CORRELATION_SAMPLE_METADATA_FILE="sample_metadata.csv"
export CORRELATION_DEMOGRAPHICS_FILE="sample_demographics.xlsx"
export CORRELATION_PROTEIN_METADATA_FILE="protein_metadata.csv"
export CORRELATION_TAU181_FILE="plasma_tau181_abeta.xlsx"
export CORRELATION_TAU217_FILE="plasma_tau217_nfl_gfap.xlsx"
export MANUSCRIPT_PANEL_SOMAKEYS=/path/to/manuscript_panel_somakeys.csv
```

If the biomarker spreadsheet needs a local study/cohort filter, also set `BIOMARKER_STUDY_FILTER`.

For pathway analyses, provide manuscript gene-list text files with one HGNC symbol per line and set:

``` bash
export MANUSCRIPT_BONE_GENE_LIST=/path/to/bone_gene_list.txt
export MANUSCRIPT_SENESCENCE_GENE_LIST=/path/to/senescence_gene_list.txt
```

## Run Order

Run proteome association models first after copying `proteome_association/config_template.yml` to local ignored `proteome_association/config.yml`:

``` bash
Rscript proteome_association/01_run_proteome_association.R --config=proteome_association/config.yml --group_mode=combined --sex_filter=all
Rscript proteome_association/01_run_proteome_association.R --config=proteome_association/config.yml --group_mode=combined --sex_filter=F
Rscript proteome_association/01_run_proteome_association.R --config=proteome_association/config.yml --group_mode=combined --sex_filter=M
Rscript proteome_association/01_run_proteome_association.R --config=proteome_association/config.yml --group_mode=individual --sex_filter=all
Rscript proteome_association/01_run_proteome_association.R --config=proteome_association/config.yml --group_mode=individual --sex_filter=F
Rscript proteome_association/01_run_proteome_association.R --config=proteome_association/config.yml --group_mode=individual --sex_filter=M
```

Then generate the final correlation heatmaps:

``` bash
bash correlation/run_final_heatmaps_from_existing_outputs.sh
```

Run pathway enrichment after the proteome association outputs are available:

``` bash
Rscript pathway_analysis/gsea_pipeline.R
```

Install required R packages in your local environment as needed; no package installation script is included.

## Key R Packages

`tidyverse`, `readr`, `readxl`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`, `ggplot2`, `ggrepel`, `lmerTest`, `broom`, `broom.mixed`, `variancePartition`, `BiocParallel`, `limma`, `ComplexHeatmap`, `circlize`, `pheatmap`, `fgsea`, `msigdbr`, `patchwork`, `janitor`, `yaml`, `lubridate`, `forcats`, `scales`, and `knitr`.
