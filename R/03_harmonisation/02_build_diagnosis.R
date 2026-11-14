# TCGA COAD/READ -> OMOP | Phase 5: harmonised Diagnosis extraction
# Run from repository root: source("R/03_harmonisation/02_build_diagnosis.R")
# Builds a source-faithful diagnosis table; this is NOT OMOP CONDITION_OCCURRENCE.
# Requires the validated Patient extraction from 01_build_patient.R.
# Source: GDC Release 46.0 clinical JSON, not the flattened clinical.tsv.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/diagnosis.csv"
validation_file <- "results/validation/02_diagnosis_validation.csv"
coverage_file <- "results/validation/02_diagnosis_field_coverage.csv"

if (!file.exists(source_file)) stop("Missing GDC JSON: ", source_file)
if (!file.exists(patient_file)) {
  stop("Missing Patient table. Run R/03_harmonisation/01_build_patient.R first.")
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")
stopifnot("patient_id" %in% names(patient))
if (nrow(patient) != 633L || anyNA(patient$patient_id) ||
    anyDuplicated(patient$patient_id)) {
  stop("Patient table must have 633 rows with distinct non-missing patient_id values.")
}

# STEP 1. Read one record per diagnosis and remember its source JSON location.
# The pointer is relative to the CASE object. Its array index is zero-based.
records <- list()
for (i in seq_along(cases)) {
  this_case <- cases[[i]]
  patient_id <- this_case$case_id
  diagnoses <- this_case$diagnoses
  if (is.null(diagnoses)) next
  if (!is.list(diagnoses)) stop("Invalid diagnoses array for case index ", i)

  for (j in seq_along(diagnoses)) {
    d <- diagnoses[[j]]
    if (!is.list(d) || is.null(names(d))) {
      stop("Invalid diagnosis object at case ", i, ", diagnosis ", j)
    }
    records[[length(records) + 1L]] <- list(
      patient_id = patient_id,
      source_json_pointer = paste0("/diagnoses/", j - 1L),
      diagnosis = d
    )
  }
}
if (length(records) != 878L) {
  stop("Expected 878 GDC diagnosis records; found ", length(records))
}

# STEP 2. Retain ALL direct diagnosis fields. Only treatments/pathology_details
# remain nested, because their own related tables will be built separately.
child_fields <- c("treatments", "pathology_details")
all_fields <- sort(unique(unlist(lapply(records, function(r) {
  names(r$diagnosis)
}), use.names = FALSE)))
scalar_fields <- setdiff(all_fields, child_fields)

# Fail rather than silently dropping any unexpected nested source object.
for (r in records) {
  d <- r$diagnosis
  for (k in scalar_fields) {
    if (!(k %in% names(d)) || is.null(d[[k]])) next
    if (is.list(d[[k]]) || length(d[[k]]) != 1L) {
      stop("Unexpected nested/non-scalar diagnosis field: ", k,
           ". Review source model before proceeding.")
    }
  }
  for (k in intersect(child_fields, names(d))) {
    if (!is.null(d[[k]]) && !is.list(d[[k]])) {
      stop("Expected array for nested field ", k)
    }
  }
}

source_text <- function(d, k) {
  if (!(k %in% names(d)) || is.null(d[[k]])) return(NA_character_)
  x <- d[[k]]
  if (length(x) != 1L || is.na(x)) return(NA_character_)
  if (is.logical(x)) return(tolower(as.character(x)))
  as.character(x)                 # Do not recode source values prematurely.
}

# Diagnosis_id and patient_id are operational relational keys. Direct source
# attributes receive the prefix diagnosis_ (same convention as Patient).
rows <- lapply(records, function(r) {
  d <- r$diagnosis
  c(diagnosis_id = source_text(d, "diagnosis_id"),
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(scalar_fields, function(k) source_text(d, k),
                            character(1)), paste0("diagnosis_", scalar_fields)))
})
diagnosis <- as.data.frame(do.call(rbind, rows),
                           stringsAsFactors = FALSE, check.names = FALSE)
rownames(diagnosis) <- NULL

# STEP 3. Validate grain, UUIDs and parent relationships BEFORE writing files.
if (nrow(diagnosis) != 878L || anyNA(diagnosis$diagnosis_id) ||
    any(!nzchar(diagnosis$diagnosis_id)) ||
    anyDuplicated(diagnosis$diagnosis_id)) {
  stop("Diagnosis row count or diagnosis_id uniqueness/availability failed.")
}
if (anyNA(diagnosis$patient_id) ||
    any(!(diagnosis$patient_id %in% patient$patient_id))) {
  stop("Diagnosis-to-Patient foreign key validation failed.")
}
if (!identical(sort(unique(diagnosis$patient_id)),
               sort(unique(patient$patient_id)))) {
  stop("Not all 633 source patients are represented in Diagnosis.")
}
if (anyDuplicated(paste(diagnosis$patient_id,
                        diagnosis$source_json_pointer, sep = "|"))) {
  stop("Duplicate case / diagnosis JSON pointer found.")
}

# Inventory of nested children: not materialised here, but conserved in JSON.
count_nested <- function(key) sum(vapply(records, function(r) {
  x <- r$diagnosis[[key]]
  if (is.null(x)) 0L else length(x)
}, integer(1)))
treatment_count <- count_nested("treatments")
pathology_count <- count_nested("pathology_details")
if (treatment_count != 2354L || pathology_count != 631L) {
  stop("Nested record totals differ from the verified source: treatments=2354; ",
       "pathology_details=631. Found ", treatment_count, " and ", pathology_count)
}

# STEP 4. Aggregate-only validation reports (safe to review before Git commit).
validation <- data.frame(
  metric = c("diagnosis_rows", "unique_diagnosis_ids", "missing_diagnosis_ids",
             "orphan_diagnosis_patient_links", "patients_with_diagnosis",
             "diagnosis_scalar_source_fields", "nested_treatment_records_preserved_in_json",
             "nested_pathology_records_preserved_in_json"),
  value = c(nrow(diagnosis), length(unique(diagnosis$diagnosis_id)),
            sum(is.na(diagnosis$diagnosis_id)),
            sum(!(diagnosis$patient_id %in% patient$patient_id)),
            length(unique(diagnosis$patient_id)), length(scalar_fields),
            treatment_count, pathology_count),
  stringsAsFactors = FALSE
)

# Field-level presence inventory, without recording individual patient values.
coverage <- do.call(rbind, lapply(scalar_fields, function(k) {
  present <- vapply(records, function(r) k %in% names(r$diagnosis), logical(1))
  explicit_null <- vapply(records, function(r) {
    k %in% names(r$diagnosis) && is.null(r$diagnosis[[k]])
  }, logical(1))
  data.frame(source_field = paste0("root[].diagnoses[].", k),
             output_column = paste0("diagnosis_", k),
             present_records = sum(present),
             absent_records = sum(!present),
             explicit_null_records = sum(explicit_null),
             stringsAsFactors = FALSE)
}))

# STEP 5. Write only after passing checks. Patient-level output stays under data/.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(diagnosis, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 diagnosis extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Nested treatments and pathology records are NOT yet separate tables.\n")
cat("This is NOT an OMOP CONDITION_OCCURRENCE table.\n")
