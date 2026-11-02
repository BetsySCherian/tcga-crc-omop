# ==============================================================================
# 05_source_variable_comparison.R
# Purpose : Reconciliation step 3. Compare shared clinical variables between
#           GDC Release 46.0, TCGA-CDR and cBioPortal PanCanAtlas, only where
#           field definitions are compatible.
#           Each comparison is classified as: agree / disagree / missing in one
#           source / missing in both. Definitional differences (units, coding,
#           granularity) are handled explicitly before comparing.
# Note    : GDC values are taken from the diagnosis labelled 'primary' FOR
#           COMPARISON ONLY. Index-diagnosis selection for the ETL is a
#           separate, future decision.
# Output  : aggregate tables in results/exploration/ (prefix 05_)
# ==============================================================================

library(jsonlite)
library(readxl)
library(readr)
suppressPackageStartupMessages(library(dplyr))
library(tidyr)

paths <- list(
  gdc_json = "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json",
  cdr      = "data/raw/cdr/TCGA-CDR-SupplementalTableS1.xlsx",
  cbio_pat = "data/raw/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
)
out_dir <- "results/exploration"

# ------------------------------------------------------------------------------
# 0. Helpers
# ------------------------------------------------------------------------------
# Source-specific "no value" markers. They are counted first (section 2), then
# converted to NA so they can never be mistaken for real values.
missing_markers <- c("", "NA", "NAN", "NONE", "'--", "NOT REPORTED", "UNKNOWN",
                     "[NOT AVAILABLE]", "[UNKNOWN]", "[NOT APPLICABLE]",
                     "[NOT EVALUATED]", "[DISCREPANCY]")

clean_text <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | x %in% missing_markers] <- NA
  x
}
clean_num <- function(x) suppressWarnings(as.numeric(clean_text(x)))

read_cbio <- function(path) {
  n_meta <- sum(startsWith(readLines(path, n = 10), "#"))
  read_tsv(path, skip = n_meta, col_types = cols(.default = col_character()),
           na = character(), show_col_types = FALSE)   # keep raw markers
}

# Stage group (I-IV) from a sub-stage such as "STAGE IIIB"
stage_group <- function(x) sub("^STAGE (IV|III|II|I|0)[A-C]?$", "\\1", x)

# ICD-O-3 topography codes for the GDC site labels found in this cohort
icdo3_site <- c(
  "CECUM" = "C18.0", "ASCENDING COLON" = "C18.2", "HEPATIC FLEXURE OF COLON" = "C18.3",
  "TRANSVERSE COLON" = "C18.4", "SPLENIC FLEXURE OF COLON" = "C18.5",
  "DESCENDING COLON" = "C18.6", "SIGMOID COLON" = "C18.7", "COLON, NOS" = "C18.9",
  "RECTOSIGMOID JUNCTION" = "C19.9", "RECTUM, NOS" = "C20.9",
  "CONNECTIVE, SUBCUTANEOUS AND OTHER SOFT TISSUES OF ABDOMEN" = "C49.4",
  "UNKNOWN PRIMARY SITE" = "C80.9"
)

# ------------------------------------------------------------------------------
# 1. Extract comparable fields: one row per patient per source (raw text kept)
# ------------------------------------------------------------------------------
gdc <- fromJSON(paths$gdc_json, simplifyVector = TRUE)

first_primary <- function(d, field) {
  d <- d[d$classification_of_tumor %in% "primary", , drop = FALSE]
  if (nrow(d) == 0 || !field %in% names(d)) NA_character_ else as.character(d[[field]][1])
}
max_followup <- function(f) {
  if (is.null(f) || NROW(f) == 0) return(NA_character_)
  v <- suppressWarnings(as.numeric(f$days_to_follow_up))
  if (all(is.na(v))) NA_character_ else as.character(max(v, na.rm = TRUE))
}
diag_field <- function(field) sapply(gdc$diagnoses, first_primary, field = field)

gdc_raw <- tibble(
  patient_id     = gdc$submitter_id,
  sex            = gdc$demographic$sex_at_birth,
  days_to_birth  = as.character(gdc$demographic$days_to_birth),
  age_obfuscated = as.character(gdc$demographic$age_is_obfuscated),
  age_at_index   = as.character(gdc$demographic$age_at_index),
  age_at_dx_days = diag_field("age_at_diagnosis"),
  site_origin    = diag_field("tissue_or_organ_of_origin"),
  site_resection = diag_field("site_of_resection_or_biopsy"),
  histology      = diag_field("morphology"),
  stage          = diag_field("ajcc_pathologic_stage"),
  edition        = diag_field("ajcc_staging_system_edition"),
  path_t         = diag_field("ajcc_pathologic_t"),
  path_n         = diag_field("ajcc_pathologic_n"),
  path_m         = diag_field("ajcc_pathologic_m"),
  vital_status   = gdc$demographic$vital_status,
  days_to_death  = as.character(gdc$demographic$days_to_death),
  last_fu_diag   = diag_field("days_to_last_follow_up"),
  last_fu_max    = sapply(gdc$follow_ups, max_followup)
)

cdr_raw <- read_excel(paths$cdr, sheet = "TCGA-CDR", guess_max = 20000,
                      .name_repair = "unique_quiet") |>
  filter(type %in% c("COAD", "READ")) |>
  transmute(patient_id    = bcr_patient_barcode,
            sex           = gender,
            days_to_birth = as.character(birth_days_to),
            age_years     = as.character(age_at_initial_pathologic_diagnosis),
            histology_txt = histological_type,
            stage         = ajcc_pathologic_tumor_stage,
            vital_status  = vital_status,
            days_to_death = as.character(death_days_to),
            last_fu       = as.character(last_contact_days_to),
            os_status     = as.character(OS),
            os_days       = as.character(OS.time))

cbio_raw <- read_cbio(paths$cbio_pat) |>
  transmute(patient_id    = PATIENT_ID,
            sex           = SEX,
            days_to_birth = DAYS_TO_BIRTH,
            age_years     = AGE,
            site_code     = ICD_O_3_SITE,
            histology     = ICD_O_3_HISTOLOGY,
            stage         = AJCC_PATHOLOGIC_TUMOR_STAGE,
            edition       = AJCC_STAGING_EDITION,
            path_t        = PATH_T_STAGE,
            path_n        = PATH_N_STAGE,
            path_m        = PATH_M_STAGE,
            os_status     = OS_STATUS,
            last_fu       = DAYS_LAST_FOLLOWUP,
            os_months     = OS_MONTHS)

# ------------------------------------------------------------------------------
# 2. Missingness BEFORE cleaning: which raw markers does each source use?
# ------------------------------------------------------------------------------
marker_counts <- function(df, source) {
  df |>
    mutate(across(-patient_id, as.character)) |>
    pivot_longer(-patient_id, names_to = "field", values_to = "raw") |>
    mutate(marker = toupper(trimws(raw)),
           marker = ifelse(is.na(raw), "<NA>", marker)) |>
    filter(marker == "<NA>" | marker %in% missing_markers) |>
    count(field, marker = ifelse(marker == "", "<empty>", marker), name = "n") |>
    mutate(source = source, .before = 1)
}
missingness_raw <- bind_rows(marker_counts(gdc_raw, "GDC R46"),
                             marker_counts(cdr_raw, "TCGA-CDR"),
                             marker_counts(cbio_raw, "cBioPortal"))

# ------------------------------------------------------------------------------
# 3. Harmonise definitions so that like is compared with like
# ------------------------------------------------------------------------------
gdc_h <- gdc_raw |> transmute(
  patient_id,
  sex            = clean_text(sex),
  days_to_birth  = clean_num(days_to_birth),
  age_years_365_25 = floor(clean_num(age_at_dx_days) / 365.25), # days -> years, convention A
  age_years        = floor(clean_num(age_at_dx_days) / 365),    # days -> years, convention B
  age_obfuscated   = toupper(as.character(age_obfuscated)) == "TRUE",
  age_at_index     = clean_num(age_at_index),                    # convention C: GDC's own field
  site_origin    = unname(icdo3_site[clean_text(site_origin)]),  # text -> ICD-O-3 code
  site_resection = unname(icdo3_site[clean_text(site_resection)]),
  histology      = clean_text(histology),
  mucinous       = ifelse(is.na(clean_text(histology)), NA, clean_text(histology) == "8480/3"),
  stage          = clean_text(stage),
  stage_group    = stage_group(clean_text(stage)),
  edition        = clean_text(edition),
  path_t = clean_text(path_t), path_n = clean_text(path_n), path_m = clean_text(path_m),
  vital_status   = clean_text(vital_status),
  days_to_death  = clean_num(days_to_death),
  last_fu        = clean_num(last_fu_diag),
  last_fu_max    = clean_num(last_fu_max)
)

cdr_h <- cdr_raw |> transmute(
  patient_id,
  sex           = clean_text(sex),
  days_to_birth = clean_num(days_to_birth),
  age_years     = clean_num(age_years),
  mucinous      = ifelse(is.na(clean_text(histology_txt)), NA,
                         grepl("MUCINOUS", clean_text(histology_txt))),
  stage         = clean_text(stage),
  stage_group   = stage_group(clean_text(stage)),
  vital_status  = clean_text(vital_status),
  days_to_death = clean_num(days_to_death),
  last_fu       = clean_num(last_fu),
  os_days       = clean_num(os_days)
)

cbio_h <- cbio_raw |> transmute(
  patient_id,
  sex           = clean_text(sex),
  days_to_birth = clean_num(days_to_birth),
  age_years     = clean_num(age_years),
  site_code     = clean_text(site_code),
  histology     = clean_text(histology),
  mucinous      = ifelse(is.na(clean_text(histology)), NA, clean_text(histology) == "8480/3"),
  stage         = clean_text(stage),
  stage_group   = stage_group(clean_text(stage)),
  edition       = clean_text(edition),
  path_t = clean_text(path_t), path_n = clean_text(path_n), path_m = clean_text(path_m),
  vital_status  = recode(clean_text(os_status), "0:LIVING" = "ALIVE", "1:DECEASED" = "DEAD"),
  last_fu       = clean_num(last_fu),
  os_months     = clean_num(os_months)
)

# ------------------------------------------------------------------------------
# 4. Unit check: cBioPortal OS_MONTHS vs TCGA-CDR OS.time (days)
# ------------------------------------------------------------------------------
os_units <- inner_join(cdr_h  |> select(patient_id, os_days),
                       cbio_h |> select(patient_id, os_months), by = "patient_id") |>
  filter(!is.na(os_days), !is.na(os_months), os_months > 0) |>
  mutate(days_per_month = os_days / os_months)

unit_check <- tibble(
  n_pairs                    = nrow(os_units),
  min_days_per_month         = min(os_units$days_per_month),
  median_days_per_month      = median(os_units$days_per_month),
  max_days_per_month         = max(os_units$days_per_month),
  candidate_365_25_div_12    = 365.25 / 12,
  candidate_365_div_12       = 365 / 12,
  n_within_half_day_365_div_12 = sum(abs(os_units$os_months * 365 / 12 - os_units$os_days) <= 0.5)
)
# Convert cBioPortal months back to days with the verified factor
cbio_h <- cbio_h |> mutate(os_days = os_months * 365 / 12)

# ------------------------------------------------------------------------------
# 5. Compare one variable between two sources (patients present in both only)
# ------------------------------------------------------------------------------
compare_var <- function(a, b, var_a, var_b = var_a, name_a, name_b, label, tol = 0) {
  x <- inner_join(a |> select(patient_id, va = all_of(var_a)),
                  b |> select(patient_id, vb = all_of(var_b)), by = "patient_id")
  both <- !is.na(x$va) & !is.na(x$vb)
  same <- if (is.numeric(x$va)) abs(x$va - x$vb) <= tol else x$va == x$vb
  tibble(variable = label, source_a = name_a, source_b = name_b,
         n_in_both_sources = nrow(x),
         agree             = sum(both & same, na.rm = TRUE),
         disagree          = sum(both & !same, na.rm = TRUE),
         missing_a_only    = sum(is.na(x$va) & !is.na(x$vb)),
         missing_b_only    = sum(!is.na(x$va) & is.na(x$vb)),
         missing_both      = sum(is.na(x$va) & is.na(x$vb)))
}

G <- "GDC R46"; D <- "TCGA-CDR"; B <- "cBioPortal"
comparison <- bind_rows(
  compare_var(gdc_h, cdr_h,  "sex", name_a = G, name_b = D, label = "sex"),
  compare_var(gdc_h, cbio_h, "sex", name_a = G, name_b = B, label = "sex"),
  compare_var(gdc_h, cdr_h,  "days_to_birth", name_a = G, name_b = D, label = "days_to_birth"),
  compare_var(gdc_h, cbio_h, "days_to_birth", name_a = G, name_b = B, label = "days_to_birth"),
  compare_var(gdc_h, cdr_h,  "age_years_365_25", "age_years", G, D, "age_years (GDC days/365.25)"),
  compare_var(gdc_h, cdr_h,  "age_years", name_a = G, name_b = D, label = "age_years (GDC days/365)"),
  compare_var(gdc_h, cdr_h,  "age_at_index", "age_years", G, D, "age_years (GDC age_at_index)"),
  compare_var(gdc_h, cbio_h, "age_years", name_a = G, name_b = B, label = "age_years (GDC days/365)"),
  compare_var(cdr_h, cbio_h, "age_years", name_a = D, name_b = B, label = "age_years"),
  compare_var(gdc_h, cbio_h, "site_origin",    "site_code", G, B, "site ICD-O-3 (GDC tissue_or_organ_of_origin)"),
  compare_var(gdc_h, cbio_h, "site_resection", "site_code", G, B, "site ICD-O-3 (GDC site_of_resection_or_biopsy)"),
  compare_var(gdc_h, cbio_h, "histology", name_a = G, name_b = B, label = "histology ICD-O-3"),
  compare_var(gdc_h, cdr_h,  "mucinous",  name_a = G, name_b = D, label = "mucinous (8480/3 vs CDR text)"),
  compare_var(gdc_h, cdr_h,  "stage", name_a = G, name_b = D, label = "AJCC pathologic sub-stage"),
  compare_var(gdc_h, cbio_h, "stage", name_a = G, name_b = B, label = "AJCC pathologic sub-stage"),
  compare_var(cdr_h, cbio_h, "stage", name_a = D, name_b = B, label = "AJCC pathologic sub-stage"),
  compare_var(gdc_h, cdr_h,  "stage_group", name_a = G, name_b = D, label = "AJCC stage group (I-IV)"),
  compare_var(gdc_h, cbio_h, "stage_group", name_a = G, name_b = B, label = "AJCC stage group (I-IV)"),
  compare_var(gdc_h, cbio_h, "edition", name_a = G, name_b = B, label = "AJCC edition"),
  compare_var(gdc_h, cbio_h, "path_t", name_a = G, name_b = B, label = "pathologic T"),
  compare_var(gdc_h, cbio_h, "path_n", name_a = G, name_b = B, label = "pathologic N"),
  compare_var(gdc_h, cbio_h, "path_m", name_a = G, name_b = B, label = "pathologic M"),
  compare_var(gdc_h, cdr_h,  "vital_status", name_a = G, name_b = D, label = "vital status"),
  compare_var(gdc_h, cbio_h, "vital_status", name_a = G, name_b = B, label = "vital status"),
  compare_var(cdr_h, cbio_h, "vital_status", name_a = D, name_b = B, label = "vital status"),
  compare_var(gdc_h, cdr_h,  "days_to_death", name_a = G, name_b = D, label = "days to death"),
  compare_var(gdc_h, cdr_h,  "last_fu", name_a = G, name_b = D, label = "days to last follow-up (GDC diagnosis field)"),
  compare_var(gdc_h, cdr_h,  "last_fu_max", "last_fu", G, D, "days to last follow-up (GDC max follow_up record)"),
  compare_var(cdr_h, cbio_h, "last_fu", name_a = D, name_b = B, label = "days to last follow-up"),
  compare_var(cdr_h, cbio_h, "os_days", name_a = D, name_b = B, label = "OS time in days (cBio months x 365/12)", tol = 0.5)
)

# ------------------------------------------------------------------------------
# 6. What do the disagreements look like? (aggregate value pairs, no patient IDs)
# ------------------------------------------------------------------------------
disagreement_pairs <- function(a, b, var_a, var_b = var_a, name_a, name_b, label) {
  inner_join(a |> select(patient_id, value_a = all_of(var_a)),
             b |> select(patient_id, value_b = all_of(var_b)), by = "patient_id") |>
    filter(!is.na(value_a), !is.na(value_b), value_a != value_b) |>
    mutate(across(c(value_a, value_b), as.character)) |>
    count(value_a, value_b, name = "n_patients", sort = TRUE) |>
    mutate(variable = label, source_a = name_a, source_b = name_b, .before = 1)
}
disagreements <- bind_rows(
  disagreement_pairs(gdc_h, cdr_h,  "stage", name_a = G, name_b = D, label = "AJCC pathologic sub-stage"),
  disagreement_pairs(gdc_h, cbio_h, "stage", name_a = G, name_b = B, label = "AJCC pathologic sub-stage"),
  disagreement_pairs(gdc_h, cbio_h, "edition", name_a = G, name_b = B, label = "AJCC edition"),
  disagreement_pairs(gdc_h, cbio_h, "site_origin", "site_code", G, B, "site ICD-O-3 (tissue_or_organ_of_origin)"),
  disagreement_pairs(gdc_h, cbio_h, "histology", name_a = G, name_b = B, label = "histology ICD-O-3"),
  disagreement_pairs(gdc_h, cdr_h,  "vital_status", name_a = G, name_b = D, label = "vital status"),
  disagreement_pairs(cdr_h, cbio_h, "vital_status", name_a = D, name_b = B, label = "vital status"),
  disagreement_pairs(gdc_h, cdr_h,  "sex", name_a = G, name_b = D, label = "sex"),
  disagreement_pairs(gdc_h, cbio_h, "path_t", name_a = G, name_b = B, label = "pathologic T"),
  disagreement_pairs(gdc_h, cbio_h, "path_n", name_a = G, name_b = B, label = "pathologic N"),
  disagreement_pairs(gdc_h, cbio_h, "path_m", name_a = G, name_b = B, label = "pathologic M")
)
# Numeric disagreements: summarise the size of the differences
numeric_diffs <- function(a, b, var_a, var_b = var_a, name_a, name_b, label) {
  inner_join(a |> select(patient_id, va = all_of(var_a)),
             b |> select(patient_id, vb = all_of(var_b)), by = "patient_id") |>
    filter(!is.na(va), !is.na(vb), va != vb) |>
    summarise(variable = label, source_a = name_a, source_b = name_b,
              n_disagree = n(), median_diff_a_minus_b = median(va - vb),
              min_diff = min(va - vb), max_diff = max(va - vb),
              n_within_1 = sum(abs(va - vb) <= 1))
}
numeric_disagreements <- bind_rows(
  numeric_diffs(gdc_h, cdr_h,  "age_years_365_25", "age_years", G, D, "age_years (GDC days/365.25)"),
  numeric_diffs(gdc_h, cdr_h,  "age_years", name_a = G, name_b = D, label = "age_years (GDC days/365)"),
  numeric_diffs(gdc_h, cdr_h,  "age_at_index", "age_years", G, D, "age_years (GDC age_at_index)"),
  numeric_diffs(gdc_h, cdr_h,  "days_to_birth", name_a = G, name_b = D, label = "days_to_birth"),
  numeric_diffs(gdc_h, cdr_h,  "days_to_death", name_a = G, name_b = D, label = "days_to_death"),
  numeric_diffs(gdc_h, cdr_h,  "last_fu", name_a = G, name_b = D, label = "last follow-up (GDC diagnosis field)"),
  numeric_diffs(gdc_h, cdr_h,  "last_fu_max", "last_fu", G, D, "last follow-up (GDC max follow_up)"),
  numeric_diffs(cdr_h, cbio_h, "last_fu", name_a = D, name_b = B, label = "last follow-up")
)

# Histology: CDR free text vs GDC ICD-O-3 code (definitional cross-tab)
histology_crosstab <- inner_join(gdc_h |> select(patient_id, gdc_icdo3 = histology),
                                 cdr_raw |> select(patient_id, cdr_text = histology_txt),
                                 by = "patient_id") |>
  count(cdr_text, gdc_icdo3, name = "n_patients", sort = TRUE)

# Are the 1-day birth differences confined to age-obfuscated (90+) patients?
birth_vs_obfuscation <- inner_join(gdc_h |> select(patient_id, gdc_birth = days_to_birth, age_obfuscated),
                                   cdr_h |> select(patient_id, cdr_birth = days_to_birth),
                                   by = "patient_id") |>
  filter(!is.na(gdc_birth), !is.na(cdr_birth)) |>
  count(age_obfuscated, birth_differs = gdc_birth != cdr_birth, name = "n_patients")

# ------------------------------------------------------------------------------
# 7. Patients with genuine contradictions (IDs + issue only; no clinical values)
# ------------------------------------------------------------------------------
gc <- inner_join(gdc_h, cdr_h, by = "patient_id", suffix = c("_gdc", "_cdr"))
flagged <- bind_rows(
  gc |> filter(vital_status_gdc != vital_status_cdr) |>
    transmute(patient_id, issue = "vital status: GDC DEAD vs TCGA-CDR ALIVE"),
  gc |> filter(!is.na(days_to_death_gdc), !is.na(days_to_death_cdr),
               days_to_death_gdc != days_to_death_cdr) |>
    transmute(patient_id, issue = "days to death differs (GDC vs TCGA-CDR)"),
  gc |> filter(!is.na(days_to_death_gdc), is.na(days_to_death_cdr)) |>
    transmute(patient_id, issue = "death days in GDC only")
) |> arrange(issue, patient_id)

# ------------------------------------------------------------------------------
# 8. Print and save
# ------------------------------------------------------------------------------
for (nm in c("unit_check", "comparison", "disagreements", "numeric_disagreements",
             "birth_vs_obfuscation", "histology_crosstab", "flagged", "missingness_raw")) {
  cat("\n=====", nm, "=====\n"); print(as.data.frame(get(nm)), right = FALSE)
}
write_csv(comparison,            file.path(out_dir, "05_variable_comparison_summary.csv"))
write_csv(unit_check,            file.path(out_dir, "05_os_unit_check.csv"))
write_csv(disagreements,         file.path(out_dir, "05_categorical_disagreements.csv"))
write_csv(numeric_disagreements, file.path(out_dir, "05_numeric_disagreements.csv"))
write_csv(histology_crosstab,    file.path(out_dir, "05_histology_cdr_vs_gdc.csv"))
write_csv(missingness_raw,       file.path(out_dir, "05_missing_value_markers.csv"))
write_csv(birth_vs_obfuscation,  file.path(out_dir, "05_birth_vs_age_obfuscation.csv"))
write_csv(flagged,               file.path(out_dir, "05_flagged_contradictions.csv"))
