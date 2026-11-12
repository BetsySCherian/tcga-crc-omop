# TCGA COAD/READ -> OMOP | Phase 5: harmonised Exposure extraction
# Run from the repository root:
#   source("R/03_harmonisation/06_build_exposure.R")
#
# Source: GDC Release 46.0 clinical JSON (NOT exposure.tsv).
# Grain: one row per case-level exposures[] source record.
# Relationship: exposure.patient_id -> patient.patient_id.
# This is a SOURCE-FAITHFUL intermediate table, NOT an OMOP table.
# An absent exposures[] array does NOT mean the patient has no exposure.
# An alcohol_history value of "Not Reported" is NOT a negative finding.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/exposure.csv"
validation_file <- "results/validation/06_exposure_validation.csv"
coverage_file <- "results/validation/06_exposure_field_coverage.csv"

for (f in c(source_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")

if (!is.list(cases) || length(cases) != 633L || nrow(patient) != 633L) {
  stop("Expected 633 GDC cases and 633 Patient rows. Check prior steps.")
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
  stop("Patient identifiers do not match GDC case identifiers.")
}

# STEP 1. Extract every source exposure. Keep its case UUID and JSON position.
records <- list()
for (i in seq_along(cases)) {
  exp_array <- cases[[i]]$exposures
  if (is.null(exp_array)) next
  if (!is.list(exp_array) || !is.null(names(exp_array))) {
    stop("Expected exposures[] array at source case position ", i)
  }
  for (j in seq_along(exp_array)) {
    obj <- exp_array[[j]]
    if (!is.list(obj) || is.null(names(obj))) {
      stop("Invalid exposure object at source case ", i, ", position ", j)
    }
    records[[length(records) + 1L]] <- list(
      patient_id = cases[[i]]$case_id,
      source_json_pointer = paste0("/exposures/", j - 1L),
      exposure = obj
    )
  }
}
if (length(records) != 2L) {
  stop("Expected 2 exposure source records; found ", length(records))
}

# STEP 2. Discover and keep all source fields, including administrative fields.
all_fields <- sort(unique(unlist(lapply(records, function(r) {
  names(r$exposure)
}), use.names = FALSE)))
for (r in records) {
  for (field in names(r$exposure)) {
    val <- r$exposure[[field]]
    if (!is.null(val) && (is.list(val) || length(val) != 1L)) {
      stop("Unexpected nested/non-scalar Exposure field: ", field)
    }
  }
}

source_text <- function(obj, field) {
  if (!(field %in% names(obj)) || is.null(obj[[field]])) return(NA_character_)
  val <- obj[[field]]
  if (length(val) != 1L || is.na(val)) return(NA_character_)
  if (is.logical(val)) return(tolower(as.character(val)))
  as.character(val)
}

rows <- lapply(records, function(r) {
  c(exposure_id = source_text(r$exposure, "exposure_id"),
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(all_fields, function(field) {
      source_text(r$exposure, field)
    }, character(1)), paste0("exposure__", all_fields)))
})
exposure <- as.data.frame(do.call(rbind, rows),
                          stringsAsFactors = FALSE, check.names = FALSE)
rownames(exposure) <- NULL

# STEP 3. Validate primary keys, foreign keys, and the precise JSON locations.
if (nrow(exposure) != 2L || anyNA(exposure$exposure_id) ||
    any(!nzchar(exposure$exposure_id)) || anyDuplicated(exposure$exposure_id)) {
  stop("Unexpected Exposure row count or missing/duplicate exposure_id.")
}
if (anyNA(exposure$patient_id) ||
    any(!(exposure$patient_id %in% patient$patient_id))) {
  stop("Orphan Exposure -> Patient links detected.")
}
if (anyDuplicated(paste(exposure$patient_id,
                        exposure$source_json_pointer, sep = "|"))) {
  stop("Duplicate Patient + JSON-pointer relationship detected.")
}
case_index <- match(exposure$patient_id, source_patient_ids)
source_index <- as.integer(sub("^/exposures/", "",
                               exposure$source_json_pointer)) + 1L
pointer_match <- vapply(seq_len(nrow(exposure)), function(i) {
  if (is.na(case_index[i]) || is.na(source_index[i])) return(FALSE)
  source_obj <- cases[[case_index[i]]]$exposures[[source_index[i]]]
  identical(source_obj$exposure_id, exposure$exposure_id[i])
}, logical(1))
if (!all(pointer_match)) stop("Exposure JSON-pointer validation failed.")

# STEP 4. Count record availability; do not equate absence with "no exposure".
patients_with_exposure <- length(unique(exposure$patient_id))
patients_without_exposure <- nrow(patient) - patients_with_exposure
alcohol_history <- vapply(records, function(r) {
  source_text(r$exposure, "alcohol_history")
}, character(1))
not_reported_count <- sum(alcohol_history == "Not Reported", na.rm = TRUE)
if (patients_with_exposure != 2L || patients_without_exposure != 631L ||
    not_reported_count != 2L) {
  stop("Exposure coverage or alcohol_history values differ from baseline.")
}

# Field inventory preserves the distinction: field absent vs explicit JSON null.
coverage <- do.call(rbind, lapply(all_fields, function(field) {
  present <- vapply(records, function(r) {
    field %in% names(r$exposure)
  }, logical(1))
  explicit_null <- vapply(records, function(r) {
    field %in% names(r$exposure) && is.null(r$exposure[[field]])
  }, logical(1))
  data.frame(
    source_path = paste0("root[].exposures[].", field),
    output_column = paste0("exposure__", field),
    source_shape = "scalar",
    present_records = sum(present),
    absent_records = sum(!present),
    explicit_null_records = sum(explicit_null),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "exposure_rows", "unique_exposure_ids", "missing_exposure_ids",
    "orphan_exposure_patient_links", "exposure_json_pointer_mismatches",
    "patients_with_exposure_records", "patients_without_exposure_records",
    "exposure_direct_source_fields", "alcohol_history_not_reported"
  ),
  value = c(
    nrow(exposure), length(unique(exposure$exposure_id)),
    sum(is.na(exposure$exposure_id)),
    sum(!(exposure$patient_id %in% patient$patient_id)),
    sum(!pointer_match), patients_with_exposure, patients_without_exposure,
    length(all_fields), not_reported_count
  ),
  stringsAsFactors = FALSE
)

# STEP 5. Save only after validation. Patient-level outputs stay Git-ignored.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(exposure, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 Exposure extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("No exposure record does NOT imply no exposure.\n")
cat("Alcohol history = 'Not Reported' is NOT a negative clinical finding.\n")
cat("This is NOT yet an OMOP table.\n")
