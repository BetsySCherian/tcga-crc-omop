# TCGA COAD/READ -> OMOP | Phase 5: harmonised Molecular Test extraction
# Run from the repository root:
#   source("R/03_harmonisation/08_build_molecular_test.R")
#
# Source: GDC Release 46.0 clinical JSON (not RNA-seq or somatic MAF data).
# Grain: one row per follow_ups[].molecular_tests[] source record.
# Parent: molecular_test.follow_up_id -> follow_up.follow_up_id.
# Retains patient_id and a zero-based, case-relative JSON pointer for lineage.
# These are clinical test records: NOT all are somatic or germline variant calls.
# This output is NOT an OMOP MEASUREMENT or OBSERVATION table.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
patient_file <- "data/processed/harmonised/patient.csv"
follow_up_file <- "data/processed/harmonised/follow_up.csv"
out_file <- "data/processed/harmonised/molecular_test.csv"
validation_file <- "results/validation/08_molecular_test_validation.csv"
coverage_file <- "results/validation/08_molecular_test_field_coverage.csv"

for (f in c(source_file, patient_file, follow_up_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
patient <- utils::read.csv(patient_file, stringsAsFactors = FALSE,
                           check.names = FALSE, colClasses = "character")
follow_up <- utils::read.csv(follow_up_file, stringsAsFactors = FALSE,
                             check.names = FALSE, colClasses = "character")

if (length(cases) != 633L || nrow(patient) != 633L || nrow(follow_up) != 3504L) {
  stop("Source, Patient or Follow-up counts differ from validated Release 46.0 baseline.")
}
for (item in list(list(patient, "patient_id"), list(follow_up, "follow_up_id"))) {
  keys <- item[[1]][[item[[2]]]]
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop("Missing or duplicate parent identifier: ", item[[2]])
  }
}
if (!all(c("patient_id", "source_json_pointer") %in% names(follow_up))) {
  stop("Follow-up must contain patient_id and source_json_pointer. Run script 04.")
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
  stop("Source case UUIDs and Patient identifiers disagree.")
}

# STEP 1: Read nested clinical molecular tests without changing source values.
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
    tests <- f$molecular_tests
    if (is.null(tests)) next
    if (!is.list(tests) || !is.null(names(tests))) {
      stop("Expected molecular_tests[] array at case ", i, ", follow-up ", j)
    }
    for (k in seq_along(tests)) {
      m <- tests[[k]]
      if (!is.list(m) || is.null(names(m))) {
        stop("Invalid molecular test object at case ", i,
             ", follow-up ", j, ", test ", k)
      }
      records[[length(records) + 1L]] <- list(
        patient_id = one_case$case_id,
        follow_up_id = f$follow_up_id,
        source_json_pointer = paste0("/follow_ups/", j - 1L,
                                     "/molecular_tests/", k - 1L),
        molecular_test = m
      )
    }
  }
}
if (length(records) != 886L) {
  stop("Expected 886 GDC clinical molecular-test records, found ", length(records))
}

# STEP 2: Capture every directly observed molecular-test source field.
# For this release these fields are all scalar. Stop on any unexpected nesting.
all_fields <- sort(unique(unlist(lapply(records, function(r) {
  names(r$molecular_test)
}), use.names = FALSE)))
if (length(all_fields) != 15L) {
  stop("Expected 15 distinct molecular-test source fields, found ", length(all_fields))
}
for (r in records) {
  for (field in names(r$molecular_test)) {
    value <- r$molecular_test[[field]]
    if (!is.null(value) && (is.list(value) || length(value) != 1L)) {
      stop("Unexpected nested or multi-valued molecular-test field: ", field)
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
  c(molecular_test_id = source_text(r$molecular_test, "molecular_test_id"),
    follow_up_id = r$follow_up_id,
    patient_id = r$patient_id,
    source_json_pointer = r$source_json_pointer,
    stats::setNames(vapply(all_fields, function(field) {
      source_text(r$molecular_test, field)
    }, character(1)), paste0("molecular_test_", all_fields)))
})
molecular_test <- as.data.frame(do.call(rbind, rows),
                                stringsAsFactors = FALSE, check.names = FALSE)
rownames(molecular_test) <- NULL

# STEP 3: Validate source IDs, the Follow-up parent, its patient and JSON pointer.
ids <- molecular_test$molecular_test_id
if (nrow(molecular_test) != 886L || anyNA(ids) || any(!nzchar(ids)) ||
    anyDuplicated(ids)) {
  stop("Molecular-test row count, missing IDs or duplicate IDs failed validation.")
}
if (anyNA(molecular_test$follow_up_id) ||
    any(!(molecular_test$follow_up_id %in% follow_up$follow_up_id))) {
  stop("Molecular Test -> Follow-up foreign-key validation failed.")
}
if (anyNA(molecular_test$patient_id) ||
    any(!(molecular_test$patient_id %in% patient$patient_id))) {
  stop("Molecular Test -> Patient linkage validation failed.")
}
parent_row <- match(molecular_test$follow_up_id, follow_up$follow_up_id)
patient_mismatch <- molecular_test$patient_id != follow_up$patient_id[parent_row]
if (any(patient_mismatch)) {
  stop("Molecular Test patient differs from its parent Follow-up patient.")
}
parent_pointer <- sub("/molecular_tests/[0-9]+$", "",
                      molecular_test$source_json_pointer)
pointer_mismatch <- parent_pointer != follow_up$source_json_pointer[parent_row]
if (any(pointer_mismatch)) {
  stop("Molecular Test JSON location does not match its parent Follow-up.")
}
if (anyDuplicated(paste(molecular_test$patient_id,
                        molecular_test$source_json_pointer, sep = "|"))) {
  stop("Duplicate case + molecular-test JSON pointer detected.")
}

case_idx <- match(molecular_test$patient_id, source_case_ids)
follow_idx <- as.integer(sub(
  "^/follow_ups/([0-9]+)/molecular_tests/[0-9]+$", "\\1",
  molecular_test$source_json_pointer)) + 1L
test_idx <- as.integer(sub("^.*/molecular_tests/([0-9]+)$", "\\1",
                           molecular_test$source_json_pointer)) + 1L
json_pointer_matches <- vapply(seq_len(nrow(molecular_test)), function(i) {
  if (anyNA(c(case_idx[i], follow_idx[i], test_idx[i]))) return(FALSE)
  original <- cases[[case_idx[i]]]$follow_ups[[follow_idx[i]]]$molecular_tests[[test_idx[i]]]
  identical(original$molecular_test_id, molecular_test$molecular_test_id[i])
}, logical(1))
if (!all(json_pointer_matches)) {
  stop("Molecular Test original JSON position/UUID validation failed.")
}

patients_with_tests <- length(unique(molecular_test$patient_id))
follow_ups_with_tests <- length(unique(molecular_test$follow_up_id))
if (patients_with_tests != 468L || follow_ups_with_tests != 886L) {
  stop("Molecular-test parent coverage differs from verified Release 46.0 baseline.")
}

# STEP 4: Aggregate reports distinguish source-absent from explicit null.
# Both appear as NA in the relational CSV, but source JSON is preserved.
coverage <- do.call(rbind, lapply(all_fields, function(field) {
  present <- vapply(records, function(r) {
    field %in% names(r$molecular_test)
  }, logical(1))
  explicit_null <- vapply(records, function(r) {
    field %in% names(r$molecular_test) && is.null(r$molecular_test[[field]])
  }, logical(1))
  data.frame(
    source_path = paste0("root[].follow_ups[].molecular_tests[].", field),
    output_column = paste0("molecular_test_", field),
    source_shape = "scalar",
    present_records = sum(present),
    absent_records = sum(!present),
    explicit_null_records = sum(explicit_null),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "molecular_test_rows", "unique_molecular_test_ids",
    "missing_molecular_test_ids", "orphan_molecular_test_follow_up_links",
    "orphan_molecular_test_patient_links", "molecular_test_parent_patient_mismatches",
    "molecular_test_parent_pointer_mismatches", "molecular_test_json_pointer_mismatches",
    "follow_ups_with_molecular_tests", "patients_with_molecular_tests",
    "molecular_test_direct_source_fields"
  ),
  value = c(
    nrow(molecular_test), length(unique(ids)), sum(is.na(ids) | !nzchar(ids)),
    sum(!(molecular_test$follow_up_id %in% follow_up$follow_up_id)),
    sum(!(molecular_test$patient_id %in% patient$patient_id)),
    sum(patient_mismatch), sum(pointer_mismatch), sum(!json_pointer_matches),
    follow_ups_with_tests, patients_with_tests, length(all_fields)
  ),
  stringsAsFactors = FALSE
)

# STEP 5: Write only after every test passes. Detailed data stays in ignored data/.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(molecular_test, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 Molecular Test extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Molecular tests link to Follow-up, not automatically to Diagnosis.\n")
cat("A clinical molecular-test record is NOT necessarily a somatic variant.\n")
cat("This is NOT yet an OMOP MEASUREMENT or OBSERVATION table.\n")
