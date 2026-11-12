# TCGA COAD/READ -> OMOP | Phase 5 integrated harmonisation validation
# Run from the repository root:
#   source("R/03_harmonisation/15_validate_harmonised_model.R")
#
# Purpose: independently verify the 15 harmonised tables together, including
# source-to-output key/parent relationships (clinical JSON + biospecimen TSVs).
# Aggregate-only reports; patient-level data remain in ignored data/ folders.
# This script does NOT generate OMOP tables or convert data types.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Required package missing: install.packages('jsonlite')")
}

base_dir <- "data/processed/harmonised"
source_dir <- "data/raw/gdc_r46"
report_path <- "results/validation/15_harmonised_integrated_validation.csv"

# Entity name, expected grain, primary key and parent foreign keys.
spec <- list(
  patient                  = list(n = 633L,   pk = "patient_id",                  fk = character()),
  diagnosis                = list(n = 878L,   pk = "diagnosis_id",                fk = c(patient_id = "patient")),
  treatment                = list(n = 2354L,  pk = "treatment_id",                fk = c(diagnosis_id = "diagnosis")),
  follow_up                = list(n = 3504L,  pk = "follow_up_id",                fk = c(patient_id = "patient")),
  family_history           = list(n = 539L,   pk = "family_history_id",           fk = c(patient_id = "patient")),
  exposure                 = list(n = 2L,     pk = "exposure_id",                 fk = c(patient_id = "patient")),
  pathology_detail         = list(n = 631L,   pk = "pathology_detail_id",         fk = c(diagnosis_id = "diagnosis")),
  molecular_test           = list(n = 886L,   pk = "molecular_test_id",           fk = c(follow_up_id = "follow_up")),
  other_clinical_attribute = list(n = 590L,   pk = "other_clinical_attribute_id", fk = c(follow_up_id = "follow_up")),
  sample                   = list(n = 2037L,  pk = "sample_id",                   fk = c(patient_id = "patient")),
  portion                  = list(n = 3217L,  pk = "portion_id",                  fk = c(sample_id = "sample")),
  analyte                  = list(n = 4436L,  pk = "analyte_id",                  fk = c(portion_id = "portion")),
  aliquot                  = list(n = 12574L, pk = "aliquot_id",                  fk = c(analyte_id = "analyte")),
  slide                    = list(n = 1998L,  pk = "slide_id",                    fk = character()),
  slide_portion            = list(n = 2035L,  pk = c("slide_id", "portion_id"),    fk = c(slide_id = "slide", portion_id = "portion"))
)

read_text <- function(path, sep) {
  if (!file.exists(path)) stop("Missing input: ", path, ". Run scripts 01-14 first.")
  utils::read.table(
    path, header = TRUE, sep = sep, quote = '"', comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character(),
    fileEncoding = "UTF-8-BOM"
  )
}

all_tables <- lapply(names(spec), function(nm) {
  read_text(file.path(base_dir, paste0(nm, ".csv")), ",")
})
names(all_tables) <- names(spec)

invalid_key <- function(x) is.na(x) | !nzchar(trimws(x)) | trimws(x) == "'--"
checks <- list()
add_check <- function(entity, test, observed, expected) {
  passed <- isTRUE(all.equal(as.character(observed), as.character(expected),
                             check.attributes = FALSE))
  checks[[length(checks) + 1L]] <<- data.frame(
    entity = entity, check = test, observed = as.character(observed),
    expected = as.character(expected), status = if (passed) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}
require_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(names(df))) stop("Duplicate column names: ", label)
}

# Stage A: row grain, unique primary keys and foreign keys for all tables.
for (nm in names(spec)) {
  df <- all_tables[[nm]]
  s <- spec[[nm]]
  require_columns(df, c(s$pk, names(s$fk)), nm)
  add_check(nm, "row_count", nrow(df), s$n)
  missing_pk <- Reduce(`|`, lapply(s$pk, function(col) invalid_key(df[[col]])))
  add_check(nm, "missing_primary_key", sum(missing_pk), 0L)
  # A separator impossible in UUIDs makes an unambiguous composite key.
  primary <- do.call(paste, c(df[s$pk], list(sep = "\037")))
  add_check(nm, "duplicate_primary_key", sum(duplicated(primary)), 0L)
  if (length(s$fk)) for (field in names(s$fk)) {
    parent_nm <- unname(s$fk[[field]])
    parent <- all_tables[[parent_nm]]
    parent_pk <- spec[[parent_nm]]$pk
    if (length(parent_pk) != 1L) stop("Unexpected composite parent key: ", parent_nm)
    add_check(nm, paste0("orphan_fk_", field),
              sum(invalid_key(df[[field]]) |
                  !(df[[field]] %in% parent[[parent_pk]])), 0L)
  }
}

# Stage B: independently check denormalised patient links on clinical children.
for (entry in list(
  c("treatment", "diagnosis", "diagnosis_id"),
  c("pathology_detail", "diagnosis", "diagnosis_id"),
  c("molecular_test", "follow_up", "follow_up_id"),
  c("other_clinical_attribute", "follow_up", "follow_up_id")
)) {
  child_nm <- entry[1]; parent_nm <- entry[2]; link <- entry[3]
  child <- all_tables[[child_nm]]; parent <- all_tables[[parent_nm]]
  require_columns(child, c(link, "patient_id"), child_nm)
  require_columns(parent, c(link, "patient_id"), parent_nm)
  idx <- match(child[[link]], parent[[link]])
  bad <- is.na(idx) | is.na(child$patient_id) |
         child$patient_id != parent$patient_id[idx]
  add_check(child_nm, "parent_patient_mismatch", sum(bad), 0L)
}
add_check("diagnosis", "patients_represented", length(unique(all_tables$diagnosis$patient_id)), 633L)
add_check("sample", "patients_represented", length(unique(all_tables$sample$patient_id)), 633L)

# Stage C: verify that each output UUID and direct parent relationship is
# identical to its original GDC source, without comparing clinical values yet.
# This checks preservation of RECORDS and LINKS, not full value-level fidelity.
source_links <- new.env(hash = TRUE)
for (nm in names(spec)) assign(nm, character(), envir = source_links)
encode <- function(id, parent = "") {
  if (is.null(id) || length(id) != 1L || is.na(id)) id <- "<missing>"
  if (is.null(parent) || length(parent) != 1L || is.na(parent)) parent <- "<missing>"
  paste(as.character(id), as.character(parent), sep = "\037")
}
add_source <- function(nm, id, parent = "") {
  old <- get(nm, envir = source_links, inherits = FALSE)
  assign(nm, c(old, encode(id, parent)), envir = source_links)
}
items <- function(x) if (is.null(x)) list() else x
clinical_file <- file.path(source_dir, "clinical/clinical.cohort.2026-10-29.json")
if (!file.exists(clinical_file)) stop("Missing GDC clinical JSON: ", clinical_file)
cases <- jsonlite::fromJSON(clinical_file, simplifyVector = FALSE)
for (case in cases) {
  pid <- case$case_id
  add_source("patient", pid)
  for (diagnosis in items(case$diagnoses)) {
    did <- diagnosis$diagnosis_id
    add_source("diagnosis", did, pid)
    for (treatment in items(diagnosis$treatments)) {
      add_source("treatment", treatment$treatment_id, did)
    }
    for (pathology in items(diagnosis$pathology_details)) {
      add_source("pathology_detail", pathology$pathology_detail_id, did)
    }
  }
  for (follow_up in items(case$follow_ups)) {
    fid <- follow_up$follow_up_id
    add_source("follow_up", fid, pid)
    for (test in items(follow_up$molecular_tests)) {
      add_source("molecular_test", test$molecular_test_id, fid)
    }
    for (attr in items(follow_up$other_clinical_attributes)) {
      add_source("other_clinical_attribute", attr$other_clinical_attribute_id, fid)
    }
  }
  for (fh in items(case$family_histories)) {
    add_source("family_history", fh$family_history_id, pid)
  }
  for (exposure in items(case$exposures)) {
    add_source("exposure", exposure$exposure_id, pid)
  }
}

# Biospecimen TSVs are source-faithful, with one row per association in slide.tsv.
source_spec <- list(
  sample = c("sample", "samples.sample_id", "cases.case_id"),
  portion = c("portion", "portions.portion_id", "samples.sample_id"),
  analyte = c("analyte", "analytes.analyte_id", "portions.portion_id"),
  aliquot = c("aliquot", "aliquots.aliquot_id", "analytes.analyte_id"),
  slide_portion = c("slide", "slides.slide_id", "portions.portion_id")
)
raw_slide <- NULL
for (nm in names(source_spec)) {
  info <- source_spec[[nm]]
  source_tsv <- read_text(file.path(source_dir, "biospecimen", paste0(info[1], ".tsv")), "\t")
  require_columns(source_tsv, info[2:3], paste0(info[1], ".tsv"))
  # Populate key/parent fingerprints vectorially for faster large TSV checks.
  assign(nm, paste(source_tsv[[info[2]]], source_tsv[[info[3]]], sep = "\037"),
         envir = source_links)
  if (nm == "slide_portion") raw_slide <- source_tsv
}
# Unique slide entity is built from distinct slide UUIDs.
assign("slide", paste(unique(raw_slide[["slides.slide_id"]]), "", sep = "\037"),
       envir = source_links)

# Compare sorted key+parent multisets (rather than only counts or set membership).
for (nm in names(spec)) {
  df <- all_tables[[nm]]
  s <- spec[[nm]]
  observed <- if (nm == "slide_portion") {
    paste(df$slide_id, df$portion_id, sep = "\037")
  } else {
    parent_col <- if (length(s$fk)) names(s$fk)[1] else NULL
    if (is.null(parent_col)) {
      vapply(df[[s$pk]], encode, character(1), USE.NAMES = FALSE)
    } else {
      mapply(encode, df[[s$pk]], df[[parent_col]], USE.NAMES = FALSE)
    }
  }
  original <- get(nm, envir = source_links, inherits = FALSE)
  add_check(nm, "original_source_records", length(original), s$n)
  add_check(nm, "source_key_parent_multiset_match",
            identical(sort(unname(observed)), sort(unname(original))), TRUE)
}

# Stage D: slide preservation and chain-level relationships.
add_check("slide_portion", "all_source_associations_preserved",
          nrow(all_tables$slide_portion), nrow(raw_slide))
repeats <- table(all_tables$slide_portion$slide_id)
add_check("slide_portion", "slides_with_multiple_portions", sum(repeats > 1L), 35L)
add_check("slide_portion", "distinct_slides_in_associations",
          length(unique(all_tables$slide_portion$slide_id)), 1998L)

# Save aggregate-only validation report even on failures; then stop on any FAIL.
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
report <- do.call(rbind, checks)
utils::write.csv(report, report_path, row.names = FALSE)
cat("\nPhase 5 integrated harmonisation validation\n")
cat("Entities checked:", length(spec), "\n")
cat("Checks passed:", sum(report$status == "PASS"), "of", nrow(report), "\n")
if (any(report$status == "FAIL")) {
  print(report[report$status == "FAIL", ], row.names = FALSE)
  stop("Integrated validation FAILED. See ", report_path)
}
cat("All checks PASSED. Aggregate report: ", report_path, "\n", sep = "")
cat("NOTE: This verifies source records, keys and parent links; it is not an OMOP CDM validation.\n")
