# TCGA COAD/READ -> OMOP | Phase 5: harmonised Follow-up extraction
# Run from repository root: source("R/03_harmonisation/04_build_follow_up.R")
# Requires validated Patient extraction (01_build_patient.R).
# Source: GDC Release 46.0 clinical JSON, not the flattened follow_up.tsv.
# Grain: one row per case-level follow_ups[] record.
# Follow-up is linked to Patient, NOT automatically to Diagnosis.
# Nested molecular_tests[] and other_clinical_attributes[] stay in the
# original source JSON for their own later extraction; do not flatten here.
# This is NOT an OMOP VISIT_OCCURRENCE/OBSERVATION table.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/follow_up.csv"
validation_file <- "results/validation/04_follow_up_validation.csv"
coverage_file <- "results/validation/04_follow_up_field_coverage.csv"

for (f in c(source_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")

if (!is.list(cases) || length(cases) != 633L || nrow(patient) != 633L) {
  stop("Expected 633 source cases and 633 Patient rows; verify previous steps.")
}
if (!("patient_id" %in% names(patient)) ||
    anyNA(patient$patient_id) || any(!nzchar(patient$patient_id)) ||
    anyDuplicated(patient$patient_id)) {
  stop("Patient table has missing or duplicate primary keys.")
}
source_patient_ids <- vapply(cases, function(x) {
  if (!is.list(x) || !is.character(x$case_id) || length(x$case_id) != 1L) {
    stop("Missing or invalid case_id in GDC source JSON.")
  }
  x$case_id
}, character(1))
if (anyDuplicated(source_patient_ids) ||
    !setequal(source_patient_ids, patient$patient_id)) {
  stop("Patient table case UUIDs do not match the GDC JSON.")
}

# STEP 1 — Extract follow-up objects from their original case positions.
# The JSON pointer is relative to the enclosing GDC case (zero-based index).
records <- list()
for (i in seq_along(cases)) {
  one_case <- cases[[i]]
  follow_ups <- one_case$follow_ups
  if (is.null(follow_ups)) next  # A case can legitimately have no follow-ups.
  if (!is.list(follow_ups) || !is.null(names(follow_ups))) {
    stop("Expected follow_ups[] array in case ", i)
  }
  for (j in seq_along(follow_ups)) {
    f <- follow_ups[[j]]
    if (!is.list(f) || is.null(names(f))) {
      stop("Invalid follow-up object in case ", i, ", index ", j)
    }
    records[[length(records) + 1L]] <- list(
      patient_id = one_case$case_id,
      source_json_pointer = paste0("/follow_ups/", j - 1L),
      follow_up = f
    )
  }
}
if (length(records) != 3504L) {
  stop("Expected 3,504 source follow-up records; found ", length(records))
}

# STEP 2 — Retain all direct scalar source attributes without recoding them.
# Keep nested arrays at their original JSON paths for future child tables.
child_fields <- c("molecular_tests", "other_clinical_attributes")
all_fields <- sort(unique(unlist(lapply(records, function(r) names(r$follow_up)),
                                 use.names = FALSE)))
scalar_fields <- setdiff(all_fields, child_fields)

for (r in records) {
  f <- r$follow_up
  for (key in intersect(scalar_fields, names(f))) {
    v <- f[[key]]
    if (!is.null(v) && (is.list(v) || length(v) != 1L)) {
      stop("Unexpected object/array in follow-up field: ", key,
           ". Review source before changing schema.")
    }
  }
  for (key in intersect(child_fields, names(f))) {
    v <- f[[key]]
    if (!is.null(v) && (!is.list(v) || !is.null(names(v)))) {
      stop("Expected nested array in follow-up field: ", key)
    }
  }
}

source_text <- function(obj, key) {
  if (!(key %in% names(obj)) || is.null(obj[[key]])) return(NA_character_)
  v <- obj[[key]]
  if (length(v) != 1L || is.na(v)) return(NA_character_)
  if (is.logical(v)) return(tolower(as.character(v)))
  as.character(v)
}

rows <- lapply(records, function(r) {
  f <- r$follow_up
  c(follow_up_id = source_text(f, "follow_up_id"),
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(scalar_fields,
                           function(key) source_text(f, key), character(1)),
                    paste0("follow_up__", scalar_fields)))
})
follow_up <- as.data.frame(do.call(rbind, rows),
                           stringsAsFactors = FALSE, check.names = FALSE)
rownames(follow_up) <- NULL

# STEP 3 — Validate primary key, foreign key and original JSON position.
if (nrow(follow_up) != 3504L ||
    anyNA(follow_up$follow_up_id) ||
    any(!nzchar(follow_up$follow_up_id)) ||
    anyDuplicated(follow_up$follow_up_id)) {
  stop("Missing/duplicate follow-up UUID or unexpected row count.")
}
if (anyNA(follow_up$patient_id) ||
    any(!(follow_up$patient_id %in% patient$patient_id))) {
  stop("Orphan Follow-up -> Patient foreign key detected.")
}
if (anyDuplicated(paste(follow_up$patient_id,
                        follow_up$source_json_pointer, sep = "|"))) {
  stop("Duplicate case + follow-up JSON pointer detected.")
}

# Re-read the pointer as a position within each parent case to confirm linkage.
case_index <- match(follow_up$patient_id, source_patient_ids)
source_index <- as.integer(sub("^/follow_ups/", "",
                               follow_up$source_json_pointer)) + 1L
pointer_match <- vapply(seq_len(nrow(follow_up)), function(i) {
  i_case <- case_index[i]
  i_fup <- source_index[i]
  if (is.na(i_case) || is.na(i_fup)) return(FALSE)
  source_record <- cases[[i_case]]$follow_ups[[i_fup]]
  identical(source_record$follow_up_id, follow_up$follow_up_id[i])
}, logical(1))
if (!all(pointer_match)) stop("Follow-up source JSON pointer validation failed.")

# The nested children are kept in the original source, not duplicated here.
nested_count <- function(field) sum(vapply(records, function(r) {
  child <- r$follow_up[[field]]
  if (is.null(child)) 0L else length(child)
}, integer(1)))

molecular_count <- nested_count("molecular_tests")
other_attribute_count <- nested_count("other_clinical_attributes")
if (molecular_count != 886L || other_attribute_count != 590L) {
  stop("Nested follow-up child counts differ from the Release 46.0 baseline: ",
       "molecular_tests=886, other_clinical_attributes=590. Found ",
       molecular_count, " and ", other_attribute_count)
}

# STEP 4 — Aggregate-only validation and field-coverage reports.
patients_with_follow_up <- length(unique(follow_up$patient_id))
patients_without_follow_up <- nrow(patient) - patients_with_follow_up
if (patients_with_follow_up != 628L || patients_without_follow_up != 5L) {
  stop("Follow-up patient coverage differs from the source baseline.")
}

coverage <- do.call(rbind, lapply(all_fields, function(key) {
  present <- vapply(records, function(r) key %in% names(r$follow_up),
                    logical(1))
  explicit_null <- vapply(records, function(r) {
    key %in% names(r$follow_up) && is.null(r$follow_up[[key]])
  }, logical(1))
  data.frame(
    source_path = paste0("root[].follow_ups[].", key),
    output_column = if (key %in% child_fields) NA_character_
                    else paste0("follow_up__", key),
    source_shape = if (key %in% child_fields) "child_array" else "scalar",
    present_records = sum(present),
    absent_records = sum(!present),
    explicit_null_records = sum(explicit_null),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "follow_up_rows", "unique_follow_up_ids", "missing_follow_up_ids",
    "orphan_follow_up_patient_links", "follow_up_json_pointer_mismatches",
    "patients_with_follow_up", "patients_without_follow_up",
    "follow_up_scalar_source_fields", "follow_up_nested_array_fields",
    "nested_molecular_test_records_preserved_in_json",
    "nested_other_clinical_attribute_records_preserved_in_json"
  ),
  value = c(
    nrow(follow_up), length(unique(follow_up$follow_up_id)),
    sum(is.na(follow_up$follow_up_id)),
    sum(!(follow_up$patient_id %in% patient$patient_id)),
    sum(!pointer_match), patients_with_follow_up, patients_without_follow_up,
    length(scalar_fields), length(intersect(child_fields, all_fields)),
    molecular_count, other_attribute_count
  ),
  stringsAsFactors = FALSE
)

# STEP 5 — Save only after all validations pass.
# Patient-level data goes under Git-ignored data/; reports are aggregate-only.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(follow_up, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 follow-up extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Follow-ups link directly to Patient, not automatically to Diagnosis.\n")
cat("Nested molecular_tests and other_clinical_attributes are not yet tables.\n")
cat("This is NOT an OMOP VISIT_OCCURRENCE or OBSERVATION table.\n")
