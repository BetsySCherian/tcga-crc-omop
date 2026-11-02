# ==============================================================================
# 03_source_patient_counts.R
# Purpose : Reconciliation step 1. Count patients in each source:
#           GDC Release 46.0 (primary), TCGA-CDR (survival reference),
#           cBioPortal PanCanAtlas (reconciliation only).
# Output  : results/exploration/source_patient_counts.csv  (aggregate counts only)
# ==============================================================================

library(jsonlite)
library(readxl)
library(readr)
suppressPackageStartupMessages(library(dplyr))

paths <- list(
  gdc_json = "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json",
  cdr      = "data/raw/cdr/TCGA-CDR-SupplementalTableS1.xlsx",
  cbio_pat = "data/raw/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  cbio_smp = "data/raw/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_sample.txt"
)

# Helper: one summary row per source -------------------------------------------
summarise_ids <- function(source, unit, ids, project = NA_character_) {
  tibble(
    source            = source,
    unit              = unit,
    project           = project,
    n_rows            = length(ids),
    n_unique_patients = n_distinct(ids),
    n_duplicate_rows  = length(ids) - n_distinct(ids),
    n_missing_id      = sum(is.na(ids) | ids == "")
  )
}

# 1. GDC (JSON: one element per case) --------------------------------------------
gdc <- fromJSON(paths$gdc_json, simplifyVector = TRUE)
gdc_ids <- tibble(patient_id = gdc$submitter_id,
                  project    = gdc$project$project_id)

# 2. TCGA-CDR (all cancers in one sheet -> keep COAD and READ) -------------------
cdr <- read_excel(paths$cdr, sheet = "TCGA-CDR", guess_max = 20000,
                  .name_repair = "unique_quiet")   # first column has no header
cdr_crc <- cdr |>
  filter(type %in% c("COAD", "READ")) |>
  transmute(patient_id = bcr_patient_barcode,
            project    = paste0("TCGA-", type))

counts <- bind_rows(
  summarise_ids("GDC R46",  "case",    gdc_ids$patient_id, "All"),
  summarise_ids("GDC R46",  "case",    gdc_ids$patient_id[gdc_ids$project == "TCGA-COAD"], "TCGA-COAD"),
  summarise_ids("GDC R46",  "case",    gdc_ids$patient_id[gdc_ids$project == "TCGA-READ"], "TCGA-READ"),
  summarise_ids("TCGA-CDR", "patient", cdr_crc$patient_id, "All"),
  summarise_ids("TCGA-CDR", "patient", cdr_crc$patient_id[cdr_crc$project == "TCGA-COAD"], "TCGA-COAD"),
  summarise_ids("TCGA-CDR", "patient", cdr_crc$patient_id[cdr_crc$project == "TCGA-READ"], "TCGA-READ")
)

# 3. cBioPortal (only if downloaded) ---------------------------------------------
# cBioPortal tables start with 4 metadata lines beginning with '#'.
read_cbio <- function(path) {
  n_meta <- sum(startsWith(readLines(path, n = 10), "#"))
  read_tsv(path, skip = n_meta, col_types = cols(.default = col_character()),
           na = c("", "NA"), show_col_types = FALSE)
}

if (file.exists(paths$cbio_pat) && file.exists(paths$cbio_smp)) {
  cbio_pat <- read_cbio(paths$cbio_pat)
  cbio_smp <- read_cbio(paths$cbio_smp)
  counts <- bind_rows(
    counts,
    summarise_ids("cBioPortal PanCanAtlas", "patient table", cbio_pat$PATIENT_ID, "All"),
    summarise_ids("cBioPortal PanCanAtlas", "patient table",
                  cbio_pat$PATIENT_ID[cbio_pat$CANCER_TYPE_ACRONYM == "COAD"], "TCGA-COAD"),
    summarise_ids("cBioPortal PanCanAtlas", "patient table",
                  cbio_pat$PATIENT_ID[cbio_pat$CANCER_TYPE_ACRONYM == "READ"], "TCGA-READ"),
    summarise_ids("cBioPortal PanCanAtlas", "sample table",  cbio_smp$PATIENT_ID, "All")
  )
} else {
  message("cBioPortal files not found: run R/01_acquisition/02_download_cbioportal_clinical.R first.")
}

print(counts, width = Inf)
write_csv(counts, "results/exploration/source_patient_counts.csv")
