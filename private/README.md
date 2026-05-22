# private/

Cohort-specific files that contain participant identifiers. **These are never committed** —
`.gitignore` excludes everything in this folder except the `.example` templates and this README.

Copy each template to its real filename and fill in your own values:

| Real file | Template | Purpose |
|---|---|---|
| `exclude_subject_ids.txt` | `exclude_subject_ids.txt.example` | One subject ID per line; these participants are dropped from every analysis. |
| `correlation_biomarker_visit_fixes.csv` | `correlation_biomarker_visit_fixes.csv.example` | Manual visit/draw-date corrections applied during biomarker harmonisation in `correlation/correlation_model_setup.Rmd`. |

Both files are optional — if absent, no exclusions and no visit fixes are applied.

### `correlation_biomarker_visit_fixes.csv` columns

- `source` — biomarker table the fix applies to (`tau181` or `tau217`).
- `id` — biomarker-table subject ID.
- `match_draw_date` — only rows with this draw date are corrected; leave blank to match all.
- `only_if_visit_missing` — `true` to apply only when the visit number is missing.
- `visit_n` — visit number to set (blank = leave unchanged).
- `draw_date` — draw date to set (blank = leave unchanged).
