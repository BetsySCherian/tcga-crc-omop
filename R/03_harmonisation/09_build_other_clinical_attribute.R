# TCGA COAD/READ -> OMOP | Phase 5: Other Clinical Attribute extraction
# Run from the repository root:
#   source("R/03_harmonisation/09_build_other_clinical_attribute.R")
#
# Source: GDC Release 46.0 clinical JSON.
# Grain: one row per follow_ups[].other_clinical_attributes[] source record.
# Parent: other_clinical_attribute.follow_up_id -> follow_up.follow_up_id.
# Store all observed direct fields; comorbidities[] and risk_factors[]
# remain JSON arrays in dedicated CSV columns rather than being collapsed.
# Do not infer a diagnosis relationship from timepoint_category.
# This is NOT yet an OMOP table.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
follow_up_file <- "data/processed/harmonised/follow_up.csv"
out_file <- "data/processed/harmonised/other_clinical_attribute.csv"
validation_file <- "results/validation/09_other_clinical_attribute_validation.csv"
coverage_file <- "results/validation/09_other_clinical_attribute_field_coverage.csv"
timepoint_file <- "results/validation/09_other_clinical_attribute_timepoints.csv"

for (f in c(source_file, patient_file, follow_up_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")
follow_up <- utils::read.csv(follow_up_file, stringsAsFactors = FALSE,
                             check.names = FALSE, colClasses = "character")

if (length(cases) != 633L || nrow(patient) != 633L || nrow(follow_up) != 3504L) {
  stop("Source/parent counts differ from the validated GDC Release 46.0 baseline.")
}
for (item in list(list(patient, "patient_id"), list(follow_up, "follow_up_id"))) {
  keys <- item[[1]][[item[[2]]]]
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop("Missing or duplicate parent identifier: ", item[[2]])
  }
}
if (!all(c("patient_id", "source_json_pointer") %in% names(follow_up))) {
  stop("Follow-up must include patient_id and source_json_pointer. Run script 04.")
}
source_case_ids <- vapply(cases, function(one_case) {
  id <- one_case$case_id
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("Invalid GDC case_id in source JSON")
  }
  id
}, character(1))
if (anyDuplicated(source_case_ids) ||
    !setequal(source_case_ids, patient$patient_id)) {
  stop("Source case UUIDs and harmonised Patient IDs disagree.")
}

# STEP 1 — Collect every nested attribute record and its exact source location.
records <- list()
for (i in seq_along(cases)) {
  one_case <- cases[[i]]
  follow_ups <- one_case$follow_ups
  if (is.null(follow_ups)) next
  if (!is.list(follow_ups) || !is.null(names(follow_ups))) {
    stop("Expected follow_ups[] array at case ", i)
  }
  for (j in seq_along(follow_ups)) {
    f <- follow_ups[[j]]
    if (!is.list(f) || is.null(names(f))) {
      stop("Invalid follow-up object at case ", i, ", follow-up ", j)
    }
    attrs <- f$other_clinical_attributes
    if (is.null(attrs)) next
    if (!is.list(attrs) || !is.null(names(attrs))) {
      stop("Expected other_clinical_attributes[] array at case ", i,
           ", follow-up ", j)
    }
    for (k in seq_along(attrs)) {
      a <- attrs[[k]]
      if (!is.list(a) || is.null(names(a))) {
        stop("Invalid clinical attribute object at case ", i,
             ", follow-up ", j, ", record ", k)
      }
      records[[length(records) + 1L]] <- list(
        patient_id = one_case$case_id,
        follow_up_id = f$follow_up_id,
        source_json_pointer = paste0("/follow_ups/", j - 1L,
                                     "/other_clinical_attributes/", k - 1L),
        attribute = a
      )
    }
  }
}
if (length(records) != 590L) {
  stop("Expected 590 other clinical attribute records; found ", length(records))
}

# STEP 2 — Preserve direct scalar fields plus the two genuine JSON arrays.
all_fields <- sort(unique(unlist(lapply(records, function(r) {
  names(r$attribute)
}), use.names = FALSE)))
array_fields <- c("comorbidities", "risk_factors")
scalar_fields <- setdiff(all_fields, array_fields)
if (length(all_fields) != 11L ||
    !all(array_fields %in% all_fields) || length(scalar_fields) != 9L) {
  stop("Other clinical attribute field inventory differs from source baseline.")
}
for (r in records) {
  for (field in names(r$attribute)) {
    v <- r$attribute[[field]]
    if (is.null(v)) next
    if (field %in% array_fields) {
      if (!is.list(v) || !is.null(names(v))) {
        stop("Expected JSON array for field: ", field)
      }
    } else if (is.list(v) || length(v) != 1L) {
      stop("Unexpected nested or multi-valued scalar field: ", field)
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

source_array_json <- function(obj, field) {
  if (!(field %in% names(obj)) || is.null(obj[[field]])) return(NA_character_)
  as.character(jsonlite::toJSON(obj[[field]], auto_unbox = FALSE,
                                null = "null", na = "null"))
}

rows <- lapply(records, function(r) {
  c(other_clinical_attribute_id = source_text(r$attribute,
                                               "other_clinical_attribute_id"),
    follow_up_id = r$follow_up_id,
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(scalar_fields, function(field) {
      source_text(r$attribute, field)
    }, character(1)), paste0("other_clinical_attribute_", scalar_fields)),
    stats::setNames(vapply(array_fields, function(field) {
      source_array_json(r$attribute, field)
    }, character(1)), paste0("other_clinical_attribute_", array_fields)))
})
other_clinical_attribute <- as.data.frame(do.call(rbind, rows),
                                           stringsAsFactors = FALSE,
                                           check.names = FALSE)
rownames(other_clinical_attribute) <- NULL

# STEP 3 — Key, parent, patient and source-position validation.
ids <- other_clinical_attribute$other_clinical_attribute_id
if (nrow(other_clinical_attribute) != 590L || anyNA(ids) ||
    any(!nzchar(ids)) || anyDuplicated(ids)) {
  stop("Other Clinical Attribute row count or identifier validation failed.")
}
if (anyNA(other_clinical_attribute$follow_up_id) ||
    any(!(other_clinical_attribute$follow_up_id %in% follow_up$follow_up_id))) {
  stop("Other Clinical Attribute -> Follow-up foreign-key validation failed.")
}
if (anyNA(other_clinical_attribute$patient_id) ||
    any(!(other_clinical_attribute$patient_id %in% patient$patient_id))) {
  stop("Other Clinical Attribute -> Patient linkage validation failed.")
}
parent_idx <- match(other_clinical_attribute$follow_up_id, follow_up$follow_up_id)
patient_mismatch <- other_clinical_attribute$patient_id != follow_up$patient_id[parent_idx]
if (any(patient_mismatch)) stop("Patient disagrees with parent Follow-up.")
parent_pointer <- sub("/other_clinical_attributes/[0-9]+$", "",
                      other_clinical_attribute$source_json_pointer)
pointer_mismatch <- parent_pointer != follow_up$source_json_pointer[parent_idx]
if (any(pointer_mismatch)) stop("JSON position disagrees with parent Follow-up.")
if (anyDuplicated(paste(other_clinical_attribute$patient_id,
                        other_clinical_attribute$source_json_pointer, sep = "|"))) {
  stop("Duplicate case + JSON pointer detected.")
}
case_idx <- match(other_clinical_attribute$patient_id, source_case_ids)
follow_idx <- as.integer(sub(
  "^/follow_ups/([0-9]+)/other_clinical_attributes/[0-9]+$", "\\1",
  other_clinical_attribute$source_json_pointer)) + 1L
attribute_idx <- as.integer(sub(
  "^.*/other_clinical_attributes/([0-9]+)$", "\\1",
  other_clinical_attribute$source_json_pointer)) + 1L
json_pointer_matches <- vapply(seq_len(nrow(other_clinical_attribute)), function(i) {
  if (anyNA(c(case_idx[i], follow_idx[i], attribute_idx[i]))) return(FALSE)
  original <- cases[[case_idx[i]]]$follow_ups[[follow_idx[i]]]$other_clinical_attributes[[attribute_idx[i]]]
  identical(original$other_clinical_attribute_id,
            other_clinical_attribute$other_clinical_attribute_id[i])
}, logical(1))
if (!all(json_pointer_matches)) {
  stop("Original JSON pointer/UUID validation failed.")
}

n_patients <- length(unique(other_clinical_attribute$patient_id))
n_followups <- length(unique(other_clinical_attribute$follow_up_id))
if (n_patients != 450L || n_followups != 590L) {
  stop("Parent coverage differs from verified Release 46.0 baseline.")
}

# STEP 4 — Report field coverage, distinguishing absent fields and explicit nulls.
coverage <- do.call(rbind, lapply(all_fields, function(field) {
  present <- vapply(records, function(r) field %in% names(r$attribute), logical(1))
  explicit_null <- vapply(records, function(r) {
    field %in% names(r$attribute) && is.null(r$attribute[[field]])
  }, logical(1))
  data.frame(
    source_path = paste0("root[].follow_ups[].other_clinical_attributes[].", field),
    output_column = paste0("other_clinical_attribute_", field),
    source_shape = if (field %in% array_fields) "array_json" else "scalar",
    present_records = sum(present),
    absent_records = sum(!present),
    explicit_null_records = sum(explicit_null),
    stringsAsFactors = FALSE
  )
}))

timepoint_col <- "other_clinical_attribute_timepoint_category"
timepoint_counts <- as.data.frame(table(
  other_clinical_attribute[[timepoint_col]], useNA = "ifany"),
  stringsAsFactors = FALSE)
names(timepoint_counts) <- c("timepoint_category", "records")

validation <- data.frame(
  metric = c(
    "other_clinical_attribute_rows", "unique_other_clinical_attribute_ids",
    "missing_other_clinical_attribute_ids", "orphan_attribute_follow_up_links",
    "orphan_attribute_patient_links", "attribute_parent_patient_mismatches",
    "attribute_parent_pointer_mismatches", "attribute_json_pointer_mismatches",
    "follow_ups_with_attributes", "patients_with_attributes",
    "attribute_direct_source_fields", "attribute_array_source_fields",
    "timepoint_initial_diagnosis", "timepoint_prior_to_diagnosis"
  ),
  value = c(
    nrow(other_clinical_attribute), length(unique(ids)),
    sum(is.na(ids) | !nzchar(ids)),
    sum(!(other_clinical_attribute$follow_up_id %in% follow_up$follow_up_id)),
    sum(!(other_clinical_attribute$patient_id %in% patient$patient_id)),
    sum(patient_mismatch), sum(pointer_mismatch), sum(!json_pointer_matches),
    n_followups, n_patients, length(all_fields), length(array_fields),
    sum(other_clinical_attribute[[timepoint_col]] == "Initial Diagnosis", na.rm = TRUE),
    sum(other_clinical_attribute[[timepoint_col]] == "Prior to Diagnosis", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

# STEP 5 — Only write after validations pass. Patient-level output stays in ignored data/.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(other_clinical_attribute, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)
utils::write.csv(timepoint_counts, timepoint_file, row.names = FALSE)

cat("\nPhase 5 Other Clinical Attribute extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("comorbidities and risk_factors retain their JSON array structure in CSV cells.\n")
cat("Timepoint category is NOT a verified diagnosis event date.\n")
cat("This is NOT yet an OMOP table.\n")
