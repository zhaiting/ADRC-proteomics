#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

Rscript correlation/make_fdr20_combined_group_heatmaps.R
Rscript correlation/make_fdr20_individual_group_heatmaps.R
