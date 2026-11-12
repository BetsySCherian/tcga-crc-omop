# TCGA COAD/READ -> OMOP | Phase 5: harmonised Treatment extraction
# Run from repository root: source("R/03_harmonisation/03_build_treatment.R")
# Requires Patient and Diagnosis tables from scripts 01 and 02.
# Source: GDC Release 46.0 clinical JSON (not the flattened TSV).
# One row = one GDC diagnoses[].treatments[] source record.
# A treatment record is NOT necessarily evidence of administered therapy.
# This is NOT an OMOP DRUG_EXPOSURE or PROCEDURE_OCCURRENCE table.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
diagnosis_file <- "data/processed/harmonised/diagnosis.csv"
out_file <- "data/processed/harmonised/treatment.csv"
validation_file <- "results/validation/03_treatment_validation.csv"
coverage_file <- "results/validation/03_treatment_field_coverage.csv"
status_file <- "results/validation/03_treatment_status_counts.csv"

for (f in c(source_file, patient_file, diagnosis_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")
diagnosis <- utils::read.csv(diagnosis_file, stringsAsFactors = FALSE,
                             check.names = FALSE, colClasses = "character")
if (length(cases) != 633L || nrow(patient) != 633L || nrow(diagnosis) != 878L) {
  stop("Source/parent-table count mismatch. Check the previous scripts and release.")
}
for (pair in list(list(patient, "patient_id"), list(diagnosis, "diagnosis_id"))) {
  v <- pair[[1]][[pair[[2]]]]
  if (is.null(v) || anyNA(v) || any(!nzchar(v)) || anyDuplicated(v)) {
    stop("Missing or duplicate primary key in parent table: ", pair[[2]])
  }
}
if (!all(c("patient_id", "source_json_pointer") %in% names(diagnosis))) {
  stop("The Diagnosis table must contain patient_id and source_json_pointer.")
}

# STEP 1: Discover each treatment at its original case/diagnosis position.
# JSON pointers are zero-based (per the earlier Diagnosis script).
records <- list()
for (i in seq_along(cases)) {
  c <- cases[[i]]
  if (is.null(c$diagnoses)) next
  if (!is.list(c$diagnoses)) stop("Invalid diagnoses array at case ", i)
  for (j in seq_along(c$diagnoses)) {
    d <- c$diagnoses[[j]]
    if (!is.list(d) || is.null(names(d))) stop("Invalid diagnosis object")
    treatments <- d$treatments
    if (is.null(treatments)) next
    if (!is.list(treatments)) stop("Invalid treatments array at case ", i, ", diagnosis ", j)
    for (k in seq_along(treatments)) {
      t <- treatments[[k]]
      if (!is.list(t) || is.null(names(t))) stop("Invalid treatment object")
      records[[length(records) + 1L]] <- list(
        patient_id = c$case_id,
        diagnosis_id = d$diagnosis_id,
        diagnosis_pointer = paste0("/diagnoses/", j - 1L),
        source_json_pointer = paste0("/diagnoses/", j - 1L,
                                     "/treatments/", k - 1L),
        treatment = t
      )
    }
  }
}
if (length(records) != 2354L) {
  stop("Expected 2,354 source treatment records; found ", length(records))
}

# STEP 2: Discover all direct source fields. This GDC release has two
# array-valued treatment attributes; retain each as a JSON array string.
# We do not explode array elements into additional treatment rows.
all_fields <- sort(unique(unlist(lapply(records, function(r) names(r$treatment)),
                                 use.names = FALSE)))
array_fields <- c("route_of_administration", "treatment_anatomic_sites")
scalar_fields <- setdiff(all_fields, array_fields)

for (r in records) {
  t <- r$treatment
  for (f in intersect(scalar_fields, names(t))) {
    v <- t[[f]]
    if (!is.null(v) && (is.list(v) || length(v) != 1L)) {
      stop("Unexpected non-scalar treatment field: ", f,
           ". Review the source before changing the schema.")
    }
  }
  for (f in intersect(array_fields, names(t))) {
    v <- t[[f]]
    if (is.null(v)) next
    if (!is.list(v) || !is.null(names(v))) {
      stop("Expected a JSON array for treatment field: ", f)
    }
    if (any(!vapply(v, function(z) is.character(z) && length(z) == 1L,
                    logical(1)))) {
      stop("Non-string array element in field ", f,
           "; do not silently flatten it.")
    }
  }
}

source_text <- function(obj, f) {
  if (!(f %in% names(obj)) || is.null(obj[[f]])) return(NA_character_)
  v <- obj[[f]]
  if (length(v) != 1L || is.na(v)) return(NA_character_)
  if (is.logical(v)) return(tolower(as.character(v)))
  as.character(v)
}
array_json <- function(obj, f) {
  if (!(f %in% names(obj)) || is.null(obj[[f]])) return(NA_character_)
  v <- obj[[f]]
  # auto_unbox=FALSE guarantees that even one-element values remain arrays.
  as.character(jsonlite::toJSON(unlist(v, use.names = FALSE),
                                auto_unbox = FALSE, null = "null"))
}

rows <- lapply(records, function(r) {
  t <- r$treatment
  c(treatment_id = source_text(t, "treatment_id"),
    diagnosis_id = r$diagnosis_id,
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(scalar_fields,
                           function(f) source_text(t, f), character(1)),
                    paste0("treatment__", scalar_fields)),
    stats::setNames(vapply(array_fields,
                           function(f) array_json(t, f), character(1)),
                    paste0("treatment__", array_fields, "_json")))
})
treatment <- as.data.frame(do.call(rbind, rows),
                           stringsAsFactors = FALSE, check.names = FALSE)
rownames(treatment) <- NULL

# STEP 3: Validate the primary keys and *both* parent relationships.
if (nrow(treatment) != 2354L || anyNA(treatment$treatment_id) ||
    any(!nzchar(treatment$treatment_id)) ||
    anyDuplicated(treatment$treatment_id)) {
  stop("Treatment ID missing/duplicate or source record count mismatch.")
}
if (anyNA(treatment$diagnosis_id) ||
    any(!(treatment$diagnosis_id %in% diagnosis$diagnosis_id))) {
  stop("Orphan Treatment -> Diagnosis foreign key detected.")
}
if (anyNA(treatment$patient_id) ||
    any(!(treatment$patient_id %in% patient$patient_id))) {
  stop("Orphan Treatment -> Patient foreign key detected.")
}
parent_match <- match(treatment$diagnosis_id, diagnosis$diagnosis_id)
if (any(treatment$patient_id != diagnosis$patient_id[parent_match])) {
  stop("Treatment patient_id conflicts with its linked Diagnosis's patient_id.")
}
expected_diag_pointer <- sub("/treatments/[0-9]+$", "",
                             treatment$source_json_pointer)
if (any(expected_diag_pointer != diagnosis$source_json_pointer[parent_match])) {
  stop("Treatment's JSON path disagrees with its Diagnosis parent path.")
}
if (anyDuplicated(paste(treatment$patient_id,
                        treatment$source_json_pointer, sep = "|"))) {
  stop("Duplicate case / treatment JSON pointer detected.")
}

# STEP 4: Record validation without publishing individual patient values.
therapy_status <- vapply(records, function(r) {
  source_text(r$treatment, "treatment_or_therapy")
}, character(1))
expected_status <- c(yes = 1086L, no = 1190L,
                     unknown = 76L, `not reported` = 2L)
observed_status <- table(factor(therapy_status, levels = names(expected_status)))
if (anyNA(therapy_status) ||
    !identical(unname(as.integer(observed_status)),
               unname(as.integer(expected_status)))) {
  stop("Treatment status distribution differs from Release 46.0 reference counts.")
}
status_counts <- data.frame(treatment_or_therapy = names(observed_status),
                            records = as.integer(observed_status))

coverage <- do.call(rbind, lapply(all_fields, function(f) {
  is_present <- vapply(records, function(r) f %in% names(r$treatment), logical(1))
  is_null <- vapply(records, function(r) {
    f %in% names(r$treatment) && is.null(r$treatment[[f]])
  }, logical(1))
  data.frame(source_path = paste0("root[].diagnoses[].treatments[].", f),
             output_column = if (f %in% array_fields)
               paste0("treatment__", f, "_json") else paste0("treatment__", f),
             source_shape = if (f %in% array_fields) "array" else "scalar",
             present_records = sum(is_present),
             absent_records = sum(!is_present),
             explicit_null_records = sum(is_null),
             stringsAsFactors = FALSE)
}))
validation <- data.frame(
  metric = c("treatment_rows", "unique_treatment_ids", "missing_treatment_ids",
             "orphan_treatment_diagnosis_links", "orphan_treatment_patient_links",
             "treatment_parent_patient_mismatches", "treatment_json_pointer_mismatches",
             "diagnoses_with_treatment", "patients_with_treatment",
             "treatment_direct_source_fields", "treatment_array_source_fields",
             "treatment_status_yes", "treatment_status_no",
             "treatment_status_unknown", "treatment_status_not_reported"),
  value = c(nrow(treatment), length(unique(treatment$treatment_id)),
            sum(is.na(treatment$treatment_id)),
            sum(!(treatment$diagnosis_id %in% diagnosis$diagnosis_id)),
            sum(!(treatment$patient_id %in% patient$patient_id)),
            sum(treatment$patient_id != diagnosis$patient_id[parent_match]),
            sum(expected_diag_pointer != diagnosis$source_json_pointer[parent_match]),
            length(unique(treatment$diagnosis_id)),
            length(unique(treatment$patient_id)), length(all_fields),
            length(intersect(array_fields, all_fields)),
            unname(observed_status["yes"]), unname(observed_status["no"]),
            unname(observed_status["unknown"]),
            unname(observed_status["not reported"])),
  stringsAsFactors = FALSE
)

# STEP 5: Write patient-level output only to Git-ignored data/, and aggregate
# validation reports to results/validation/. Review these reports before Git.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(treatment, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)
utils::write.csv(status_counts, status_file, row.names = FALSE)

cat("\nPhase 5 treatment extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("A treatment source record (including 'no') is NOT proof of treatment administration.\n")
cat("This is NOT an OMOP DRUG_EXPOSURE or PROCEDURE_OCCURRENCE table.\n")
