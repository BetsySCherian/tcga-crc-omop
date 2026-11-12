# TCGA COAD/READ -> OMOP | Phase 5: harmonised Pathology Detail extraction
# Run from repository root:
#   source("R/03_harmonisation/07_build_pathology_detail.R")
#
# Source: GDC Release 46.0 clinical JSON (not pathology_detail.tsv).
# Grain: one row per diagnoses[].pathology_details[] source record.
# Parent: pathology_detail.diagnosis_id -> diagnosis.diagnosis_id.
# Source lineage also retains patient_id and zero-based JSON pointer.
# This is source-faithful intermediate data, NOT an OMOP table.
# Missing pathology_details[] does NOT establish absence of pathology findings.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
diagnosis_file <- "data/processed/harmonised/diagnosis.csv"
out_file <- "data/processed/harmonised/pathology_detail.csv"
validation_file <- "results/validation/07_pathology_detail_validation.csv"
coverage_file <- "results/validation/07_pathology_detail_field_coverage.csv"

for (f in c(source_file, patient_file, diagnosis_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")
diagnosis <- utils::read.csv(diagnosis_file, stringsAsFactors = FALSE,
                             check.names = FALSE, colClasses = "character")
if (length(cases) != 633L || nrow(patient) != 633L || nrow(diagnosis) != 878L) {
  stop("Case/Patient/Diagnosis counts do not match the verified Release 46.0 baseline.")
}
for (item in list(list(patient, "patient_id"), list(diagnosis, "diagnosis_id"))) {
  keys <- item[[1]][[item[[2]]]]
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop("Missing or duplicate parent key: ", item[[2]])
  }
}
if (!all(c("patient_id", "source_json_pointer") %in% names(diagnosis))) {
  stop("Diagnosis must have patient_id and source_json_pointer. Run script 02 first.")
}

# STEP 1: Extract one record for each source pathology_detail object.
# Store its original case UUID, diagnosis UUID, and zero-based JSON location.
records <- list()
for (i in seq_along(cases)) {
  case <- cases[[i]]
  diagnoses <- case$diagnoses
  if (is.null(diagnoses)) next
  if (!is.list(diagnoses) || !is.null(names(diagnoses))) {
    stop("Invalid diagnoses[] array for case ", i)
  }
  for (j in seq_along(diagnoses)) {
    d <- diagnoses[[j]]
    if (!is.list(d) || is.null(names(d))) {
      stop("Invalid diagnosis object at case ", i, ", diagnosis ", j)
    }
    details <- d$pathology_details
    if (is.null(details)) next
    if (!is.list(details) || !is.null(names(details))) {
      stop("Invalid pathology_details[] array at case ", i, ", diagnosis ", j)
    }
    for (k in seq_along(details)) {
      p <- details[[k]]
      if (!is.list(p) || is.null(names(p))) {
        stop("Invalid pathology detail object at case ", i,
             ", diagnosis ", j, ", detail ", k)
      }
      records[[length(records) + 1L]] <- list(
        patient_id = case$case_id,
        diagnosis_id = d$diagnosis_id,
        source_json_pointer = paste0("/diagnoses/", j - 1L,
                                     "/pathology_details/", k - 1L),
        pathology_detail = p
      )
    }
  }
}
if (length(records) != 631L) {
  stop("Expected 631 pathology-detail records; found ", length(records))
}

# STEP 2: Discover every direct source attribute. Preserve original values.
# This release has scalar pathology-detail attributes only. Fail if that changes.
all_fields <- sort(unique(unlist(lapply(records, function(r) {
  names(r$pathology_detail)
}), use.names = FALSE)))
for (r in records) {
  for (field in names(r$pathology_detail)) {
    val <- r$pathology_detail[[field]]
    if (!is.null(val) && (is.list(val) || length(val) != 1L)) {
      stop("Unexpected non-scalar pathology field: ", field,
           ". Review source schema rather than flattening it.")
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
  c(pathology_detail_id = source_text(r$pathology_detail, "pathology_detail_id"),
    diagnosis_id = r$diagnosis_id,
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(all_fields, function(field) {
      source_text(r$pathology_detail, field)
    }, character(1)), paste0("pathology_detail__", all_fields)))
})
pathology_detail <- as.data.frame(do.call(rbind, rows),
                                  stringsAsFactors = FALSE, check.names = FALSE)
rownames(pathology_detail) <- NULL

# STEP 3: Validate uniqueness, parent diagnosis, patient consistency, and lineage.
ids <- pathology_detail$pathology_detail_id
if (nrow(pathology_detail) != 631L || anyNA(ids) ||
    any(!nzchar(ids)) || anyDuplicated(ids)) {
  stop("Pathology-detail row count/identifier uniqueness check failed.")
}
if (anyNA(pathology_detail$diagnosis_id) ||
    any(!(pathology_detail$diagnosis_id %in% diagnosis$diagnosis_id))) {
  stop("Pathology Detail -> Diagnosis foreign-key check failed.")
}
if (anyNA(pathology_detail$patient_id) ||
    any(!(pathology_detail$patient_id %in% patient$patient_id))) {
  stop("Pathology Detail -> Patient linkage check failed.")
}
parent_row <- match(pathology_detail$diagnosis_id, diagnosis$diagnosis_id)
patient_mismatch <- pathology_detail$patient_id != diagnosis$patient_id[parent_row]
if (any(patient_mismatch)) stop("Pathology Detail conflicts with parent Diagnosis's patient.")
parent_pointer <- sub("/pathology_details/[0-9]+$", "",
                      pathology_detail$source_json_pointer)
pointer_mismatch <- parent_pointer != diagnosis$source_json_pointer[parent_row]
if (any(pointer_mismatch)) stop("Pathology Detail JSON pointer does not match parent Diagnosis.")
if (anyDuplicated(paste(pathology_detail$patient_id,
                        pathology_detail$source_json_pointer, sep = "|"))) {
  stop("Duplicate case + pathology JSON pointer.")
}

# Verify each path resolves to the same original source ID (not only to a parent).
source_case_ids <- vapply(cases, function(c) c$case_id, character(1))
case_idx <- match(pathology_detail$patient_id, source_case_ids)
source_diag_idx <- as.integer(sub("^/diagnoses/([0-9]+)/pathology_details/[0-9]+$",
                                  "\\1", pathology_detail$source_json_pointer)) + 1L
source_detail_idx <- as.integer(sub("^.*/pathology_details/([0-9]+)$",
                                    "\\1", pathology_detail$source_json_pointer)) + 1L
json_pointer_matches <- vapply(seq_len(nrow(pathology_detail)), function(i) {
  if (anyNA(c(case_idx[i], source_diag_idx[i], source_detail_idx[i]))) return(FALSE)
  obj <- cases[[case_idx[i]]]$diagnoses[[source_diag_idx[i]]]$pathology_details[[source_detail_idx[i]]]
  identical(obj$pathology_detail_id, pathology_detail$pathology_detail_id[i])
}, logical(1))
if (!all(json_pointer_matches)) stop("Pathology-detail JSON location verification failed.")

# STEP 4: Aggregate-only completeness and validation reports.
coverage <- do.call(rbind, lapply(all_fields, function(field) {
  present <- vapply(records, function(r) {
    field %in% names(r$pathology_detail)
  }, logical(1))
  explicit_null <- vapply(records, function(r) {
    field %in% names(r$pathology_detail) && is.null(r$pathology_detail[[field]])
  }, logical(1))
  data.frame(
    source_path = paste0("root[].diagnoses[].pathology_details[].", field),
    output_column = paste0("pathology_detail__", field),
    source_shape = "scalar",
    present_records = sum(present),
    absent_records = sum(!present),
    explicit_null_records = sum(explicit_null),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "pathology_detail_rows", "unique_pathology_detail_ids",
    "missing_pathology_detail_ids", "orphan_pathology_detail_diagnosis_links",
    "orphan_pathology_detail_patient_links", "pathology_parent_patient_mismatches",
    "pathology_parent_pointer_mismatches", "pathology_json_pointer_mismatches",
    "diagnoses_with_pathology_details", "patients_with_pathology_details",
    "pathology_detail_direct_source_fields"
  ),
  value = c(
    nrow(pathology_detail), length(unique(ids)), sum(is.na(ids) | !nzchar(ids)),
    sum(!(pathology_detail$diagnosis_id %in% diagnosis$diagnosis_id)),
    sum(!(pathology_detail$patient_id %in% patient$patient_id)),
    sum(patient_mismatch), sum(pointer_mismatch), sum(!json_pointer_matches),
    length(unique(pathology_detail$diagnosis_id)),
    length(unique(pathology_detail$patient_id)), length(all_fields)
  ),
  stringsAsFactors = FALSE
)

# STEP 5: Only write after all validations pass. Detailed rows remain under Git-ignored data/.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(pathology_detail, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 Pathology Detail extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Pathology details link to Diagnosis, not independently to an inferred diagnosis event.\n")
cat("Missing pathology_details[] does NOT mean a patient had no pathology.\n")
cat("This is NOT yet an OMOP table.\n")
