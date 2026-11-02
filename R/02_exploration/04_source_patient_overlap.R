# ==============================================================================
# 04_source_patient_overlap.R
# Purpose : Reconciliation step 2. Patient-ID overlap between
#           GDC Release 46.0 (primary ETL source), TCGA-CDR (survival reference)
#           and cBioPortal PanCanAtlas (reconciliation only).
#           1) three-way and pairwise overlap   2) patients missing per source
#           3) nesting tests                     4) COAD/READ label concordance
#           5) checks of candidate explanations for discrepancies
# Output  : aggregate tables in results/exploration/ (prefix 04_)
#           + list of discrepant patients (open-access TCGA barcodes only)
# ==============================================================================

library(jsonlite)
library(readxl)
library(readr)
suppressPackageStartupMessages(library(dplyr))

paths <- list(
  gdc_json   = "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json",
  gdc_sample = "data/raw/gdc_r46/biospecimen/sample.tsv",
  cdr        = "data/raw/cdr/TCGA-CDR-SupplementalTableS1.xlsx",
  cbio_pat   = "data/raw/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
)
out_dir <- "results/exploration"

read_cbio <- function(path) {               # skip cBioPortal's '#' metadata lines
  n_meta <- sum(startsWith(readLines(path, n = 10), "#"))
  read_tsv(path, skip = n_meta, col_types = cols(.default = col_character()),
           na = c("", "NA"), show_col_types = FALSE)
}

# 1. Patient IDs and project labels from each source -----------------------------
gdc <- fromJSON(paths$gdc_json, simplifyVector = TRUE)
gdc_ids <- tibble(patient_id  = gdc$submitter_id,
                  project_gdc = gdc$project$project_id)

cdr_ids <- read_excel(paths$cdr, sheet = "TCGA-CDR", guess_max = 20000,
                      .name_repair = "unique_quiet") |>
  filter(type %in% c("COAD", "READ")) |>
  transmute(patient_id  = bcr_patient_barcode,
            project_cdr = paste0("TCGA-", type))

cbio_ids <- read_cbio(paths$cbio_pat) |>
  transmute(patient_id   = PATIENT_ID,
            project_cbio = paste0("TCGA-", CANCER_TYPE_ACRONYM))

# 2. Membership table: one row per patient found in ANY source -------------------
membership <- gdc_ids |>
  full_join(cdr_ids,  by = "patient_id") |>
  full_join(cbio_ids, by = "patient_id") |>
  mutate(in_gdc  = !is.na(project_gdc),
         in_cdr  = !is.na(project_cdr),
         in_cbio = !is.na(project_cbio))

# 3. Three-way overlap patterns ---------------------------------------------------
overlap_patterns <- membership |>
  count(in_gdc, in_cdr, in_cbio, name = "n_patients") |>
  arrange(desc(n_patients))

# 4. Pairwise overlap -------------------------------------------------------------
pairwise <- function(a, b, name_a, name_b) {
  tibble(source_a = name_a, source_b = name_b,
         n_a = sum(a), n_b = sum(b), n_both = sum(a & b),
         n_only_a = sum(a & !b), n_only_b = sum(!a & b))
}
overlap_pairwise <- with(membership, bind_rows(
  pairwise(in_gdc, in_cdr,  "GDC R46",  "TCGA-CDR"),
  pairwise(in_gdc, in_cbio, "GDC R46",  "cBioPortal"),
  pairwise(in_cdr, in_cbio, "TCGA-CDR", "cBioPortal")
))

# 5. Nesting tests: "A is within B" holds only if no patient is in A but not B ---
nesting <- with(membership, tibble(
  test         = c("TCGA-CDR within GDC", "cBioPortal within GDC", "cBioPortal within TCGA-CDR"),
  n_violations = c(sum(in_cdr & !in_gdc), sum(in_cbio & !in_gdc), sum(in_cbio & !in_cdr))
)) |> mutate(nested = n_violations == 0)

# 6. COAD/READ project-label concordance for patients present in both sources ---
project_concordance <- bind_rows(
  membership |> filter(in_gdc, in_cdr) |>
    summarise(comparison = "GDC vs TCGA-CDR", n_compared = n(),
              n_discordant = sum(project_gdc != project_cdr)),
  membership |> filter(in_gdc, in_cbio) |>
    summarise(comparison = "GDC vs cBioPortal", n_compared = n(),
              n_discordant = sum(project_gdc != project_cbio))
)

# 7. Discrepant patients: everyone not present in all three sources --------------
discrepant <- membership |>
  filter(!(in_gdc & in_cdr & in_cbio)) |>
  select(patient_id, project_gdc, in_gdc, in_cdr, in_cbio) |>
  arrange(in_cdr, in_cbio, patient_id)

# 8. Check candidate explanations -------------------------------------------------
# 8a. Do patients missing from cBioPortal lack a tumour sample in GDC?
#     (cBioPortal keeps one analysed primary-tumour sample per patient; the
#      analysed TCGA samples were frozen tissue, not FFPE blocks.)
gdc_samples <- read_tsv(paths$gdc_sample, col_types = cols(.default = col_character()),
                        na = c("", "'--"), show_col_types = FALSE)
with_nonffpe_tumour <- gdc_samples |>
  filter(samples.sample_type == "Primary Tumor",
         samples.preservation_method != "FFPE") |>
  distinct(cases.submitter_id) |> pull()

check_tumour_sample <- membership |>
  filter(in_gdc) |>
  mutate(has_nonffpe_primary_tumour_sample = patient_id %in% with_nonffpe_tumour) |>
  count(in_cbio, has_nonffpe_primary_tumour_sample, name = "n_patients")

# 8b. GDC clinical completeness for patients absent from TCGA-CDR
gdc_followups <- sapply(gdc$follow_ups, function(x) if (is.null(x)) 0L else NROW(x))
check_not_in_cdr <- gdc_ids |>
  mutate(vital_status = gdc$demographic$vital_status,
         n_follow_ups = gdc_followups) |>
  filter(!patient_id %in% cdr_ids$patient_id)

# 9. Project label vs anatomical site (GDC primary diagnosis, ICD-O-3 topography) -
primary_site <- bind_rows(lapply(seq_len(nrow(gdc)), function(i) {
  d <- gdc$diagnoses[[i]]
  d <- d[d$classification_of_tumor %in% "primary", , drop = FALSE]
  tibble(patient_id = gdc$submitter_id[i],
         site = if (nrow(d) > 0) d$tissue_or_organ_of_origin[1] else NA_character_)
}))
project_vs_site <- gdc_ids |>
  left_join(primary_site, by = "patient_id") |>
  count(project_gdc, site, name = "n_patients") |>
  arrange(project_gdc, desc(n_patients))

# 10. Print and save --------------------------------------------------------------
for (nm in c("overlap_patterns", "overlap_pairwise", "nesting", "project_concordance",
             "check_tumour_sample", "check_not_in_cdr", "project_vs_site", "discrepant")) {
  cat("\n=====", nm, "=====\n"); print(get(nm), n = Inf, width = Inf)
}
write_csv(overlap_patterns,    file.path(out_dir, "04_overlap_patterns.csv"))
write_csv(overlap_pairwise,    file.path(out_dir, "04_overlap_pairwise.csv"))
write_csv(nesting,             file.path(out_dir, "04_nesting_tests.csv"))
write_csv(project_concordance, file.path(out_dir, "04_project_concordance.csv"))
write_csv(check_tumour_sample, file.path(out_dir, "04_check_cbio_missing_tumour_sample.csv"))
write_csv(project_vs_site,     file.path(out_dir, "04_project_vs_primary_site.csv"))
write_csv(discrepant,          file.path(out_dir, "04_discrepant_patients.csv"))
