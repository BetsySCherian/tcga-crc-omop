# TCGA COAD/READ -> OMOP | Phase 5: harmonised Family History extraction
# Run from repository root:
#   source("R/03_harmonisation/05_build_family_history.R")
# Requires the validated Patient table from 01_build_patient.R.
# Source: GDC Release 46.0 clinical JSON, NOT family_history.tsv.
# Grain: one row per case-level family_histories[] source record.
# Relationship: family_history.patient_id -> patient.patient_id.
# Important: a source record with relative_with_cancer_history = "no"
# does NOT represent a positive cancer family history.
# This is NOT an OMOP FAMILY_HISTORY/OBSERVATION table.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/family_history.csv"
validation_file <- "results/validation/05_family_history_validation.csv"
coverage_file <- "results/validation/05_family_history_field_coverage.csv"
status_file <- "results/validation/05_family_history_status_counts.csv"

for (f in c(source_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")

if (!is.list(cases) || length(cases) != 633L || nrow(patient) != 633L) {
  stop("Expected 633 source cases and 633 Patient rows. Check prior steps.")
}
if (!("patient_id" %in% names(patient)) || anyNA(patient$patient_id) ||
    any(!nzchar(patient$patient_id)) || anyDuplicated(patient$patient_id)) {
  stop("Patient table has missing or duplicate primary keys.")
}
source_patient_ids <- vapply(cases, function(case) {
  id <- case$case_id
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("Missing or invalid case_id in source JSON.")
  }
  id
}, character(1))
if (anyDuplicated(source_patient_ids) ||
    !setequal(source_patient_ids, patient$patient_id)) {
  stop("Patient identifiers do not match source GDC case identifiers.")
}

# STEP 1 — Extract every family history record with its GDC case UUID
# and position within that case's family_histories[] array.
records <- list()
for (i in seq_along(cases)) {
  fam <- cases[[i]]$family_histories
  if (is.null(fam)) next
  if (!is.list(fam) || !is.null(names(fam))) {
    stop("Expected family_histories[] array at source case position ", i)
  }
  for (j in seq_along(fam)) {
    record <- fam[[j]]
    if (!is.list(record) || is.null(names(record))) {
      stop("Invalid family history object at case ", i, ", position ", j)
    }
    records[[length(records) + 1L]] <- list(
      patient_id = cases[[i]]$case_id,
      source_json_pointer = paste0("/family_histories/", j - 1L),
      family_history = record
    )
  }
}
if (length(records) != 539L) {
  stop("Expected 539 family history source records; found ", length(records))
}

# STEP 2 — Discover and preserve ALL direct source attributes.
# Source values are left as text for now, like the other Phase 5 tables.
# Field coverage separately records absent and explicit-null attributes.
all_fields <- sort(unique(unlist(lapply(records, function(r) {
  names(r$family_history)
}), use.names = FALSE)))

for (r in records) {
  for (field in names(r$family_history)) {
    value <- r$family_history[[field]]
    if (!is.null(value) && (is.list(value) || length(value) != 1L)) {
      stop("Unexpected nested/non-scalar Family History field: ", field,
           ". Inspect original JSON rather than flatten it.")
    }
  }
}

source_text <- function(obj, field) {
  if (!(field %in% names(obj)) || is.null(obj[[field]])) return(NA_character_)
  value <- obj[[field]]
  if (length(value) != 1L || is.na(value)) return(NA_character_)
  if (is.logical(value)) return(tolower(as.character(value)))
  as.character(value)
}

rows <- lapply(records, function(r) {
  record <- r$family_history
  c(family_history_id = source_text(record, "family_history_id"),
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(all_fields, function(field) {
      source_text(record, field)
    }, character(1)), paste0("family_history__", all_fields)))
})
family_history <- as.data.frame(do.call(rbind, rows),
                                stringsAsFactors = FALSE, check.names = FALSE)
rownames(family_history) <- NULL

# STEP 3 — Verify source counts, primary keys and parent-case links.
if (nrow(family_history) != 539L ||
    anyNA(family_history$family_history_id) ||
    any(!nzchar(family_history$family_history_id)) ||
    anyDuplicated(family_history$family_history_id)) {
  stop("Missing/duplicate family_history_id or unexpected row count.")
}
if (anyNA(family_history$patient_id) ||
    any(!(family_history$patient_id %in% patient$patient_id))) {
  stop("Orphan Family History -> Patient link detected.")
}
if (anyDuplicated(paste(family_history$patient_id,
                        family_history$source_json_pointer, sep = "|"))) {
  stop("Duplicate case + JSON pointer detected.")
}

# Verify each extracted identifier against its exact position in the JSON.
case_index <- match(family_history$patient_id, source_patient_ids)
source_index <- as.integer(sub("^/family_histories/", "",
                               family_history$source_json_pointer)) + 1L
pointer_match <- vapply(seq_len(nrow(family_history)), function(i) {
  if (is.na(case_index[i]) || is.na(source_index[i])) return(FALSE)
  source_obj <- cases[[case_index[i]]]$family_histories[[source_index[i]]]
  identical(source_obj$family_history_id, family_history$family_history_id[i])
}, logical(1))
if (!all(pointer_match)) stop("Family History JSON pointer validation failed.")

patients_with_family_history <- length(unique(family_history$patient_id))
patients_without_family_history <- nrow(patient) - patients_with_family_history
if (patients_with_family_history != 539L ||
    patients_without_family_history != 94L) {
  stop("Unexpected Patient coverage for Family History.")
}

# STEP 4 — Validate clinically important yes / no / unknown distinction.
family_indicator <- vapply(records, function(r) {
  source_text(r$family_history, "relative_with_cancer_history")
}, character(1))
expected_indicator <- c(yes = 74L, no = 463L, unknown = 2L)
indicator_counts <- table(factor(family_indicator,
                                 levels = names(expected_indicator)))
if (anyNA(family_indicator) ||
    !identical(unname(as.integer(indicator_counts)),
               unname(as.integer(expected_indicator)))) {
  stop("Family history cancer indicator differs from the GDC baseline.")
}
status_counts <- data.frame(
  relative_with_cancer_history = names(indicator_counts),
  records = as.integer(indicator_counts),
  stringsAsFactors = FALSE
)

coverage <- do.call(rbind, lapply(all_fields, function(field) {
  present <- vapply(records, function(r) {
    field %in% names(r$family_history)
  }, logical(1))
  explicit_null <- vapply(records, function(r) {
    field %in% names(r$family_history) &&
      is.null(r$family_history[[field]])
  }, logical(1))
  data.frame(
    source_path = paste0("root[].family_histories[].", field),
    output_column = paste0("family_history__", field),
    source_shape = "scalar",
    present_records = sum(present),
    absent_records = sum(!present),
    explicit_null_records = sum(explicit_null),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "family_history_rows", "unique_family_history_ids",
    "missing_family_history_ids", "orphan_family_history_patient_links",
    "family_history_json_pointer_mismatches",
    "patients_with_family_history_records",
    "patients_without_family_history_records",
    "family_history_direct_source_fields",
    "relative_with_cancer_history_yes",
    "relative_with_cancer_history_no",
    "relative_with_cancer_history_unknown"
  ),
  value = c(
    nrow(family_history), length(unique(family_history$family_history_id)),
    sum(is.na(family_history$family_history_id)),
    sum(!(family_history$patient_id %in% patient$patient_id)),
    sum(!pointer_match), patients_with_family_history,
    patients_without_family_history, length(all_fields),
    unname(indicator_counts["yes"]), unname(indicator_counts["no"]),
    unname(indicator_counts["unknown"])
  ),
  stringsAsFactors = FALSE
)

# STEP 5 — Only save validated outputs. data/ MUST be Git-ignored.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(family_history, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)
utils::write.csv(status_counts, status_file, row.names = FALSE)

cat("\nPhase 5 Family History extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("A 'no' record is NOT a positive cancer family history.\n")
cat("This is NOT yet an OMOP table.\n")
