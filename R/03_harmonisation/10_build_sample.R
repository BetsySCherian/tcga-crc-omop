# TCGA COAD/READ -> OMOP | Phase 5: harmonised Sample extraction
# Run from the repository root:
#   source("R/03_harmonisation/10_build_sample.R")
#
# Source: GDC Release 46.0 biospecimen/sample.tsv (NOT clinical JSON).
# Grain: one row per GDC samples.sample_id (sample UUID).
# Relationship: sample.patient_id -> patient.patient_id (GDC case UUID).
# Preserve all samples.* fields as original character values. In particular,
# the GDC TSV marker "'--" remains distinguishable from an empty cell.
# This is NOT an OMOP SPECIMEN table; no OMOP mapping is performed.

source_file <- "data/raw/gdc_r46/biospecimen/sample.tsv"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/sample.csv"
validation_file <- "results/validation/10_sample_validation.csv"
coverage_file <- "results/validation/10_sample_field_coverage.csv"
sample_types_file <- "results/validation/10_sample_types.csv"

for (f in c(source_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input: ", f, "\nRun from the project root.")
}

# STEP 1: Read the flattened GDC sample TSV and the validated Patient table.
# These files contain identifiers and source values; import as text so R does
# not silently change TCGA barcodes, numbers, leading zeros or missing markers.
raw_sample <- utils::read.delim(
  source_file, sep = "\t", header = TRUE,
  check.names = FALSE, stringsAsFactors = FALSE,
  colClasses = "character", na.strings = character(),
  quote = "\"", comment.char = "", fileEncoding = "UTF-8-BOM"
)
patient <- utils::read.csv(
  patient_file, check.names = FALSE, stringsAsFactors = FALSE,
  colClasses = "character"
)

required <- c("project.project_id", "cases.case_id", "cases.submitter_id",
              "samples.sample_id", "samples.submitter_id", "samples.sample_type")
missing_columns <- setdiff(required, names(raw_sample))
if (length(missing_columns)) {
  stop("Required biospecimen TSV columns missing: ",
       paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(names(raw_sample))) stop("Duplicate source TSV column names.")
if (!all(c("patient_id", "case_submitter_id", "project_project_id") %in%
         names(patient))) {
  stop("Expected Patient source identifiers missing. Run/check 01_build_patient.R.")
}
if (nrow(patient) != 633L || anyNA(patient$patient_id) ||
    any(!nzchar(patient$patient_id)) || anyDuplicated(patient$patient_id)) {
  stop("Patient table must have 633 distinct non-missing case UUIDs.")
}

# STEP 2: Build one source-faithful row per sample.
# The samples.* attributes (including original samples.sample_id) are retained
# with a sample_ prefix. Patient/project-level attributes belong in Patient.
source_sample_fields <- grep("^samples\\.", names(raw_sample), value = TRUE)
if (length(source_sample_fields) != 30L || ncol(raw_sample) != 33L) {
  stop("Sample TSV column inventory changed: expected 30 samples.* of 33 columns; ",
       "observed ", length(source_sample_fields), " of ", ncol(raw_sample), ".")
}

sample <- data.frame(
  sample_id = raw_sample[["samples.sample_id"]],
  patient_id = raw_sample[["cases.case_id"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (field in source_sample_fields) {
  output_name <- paste0("sample_", sub("^samples\\.", "", field))
  if (output_name %in% names(sample)) stop("Unexpected output column collision.")
  sample[[output_name]] <- raw_sample[[field]]
}

# STEP 3: Validate grain and relationships before saving.
invalid_key <- function(x) is.na(x) | !nzchar(trimws(x)) | x == "'--"
if (nrow(sample) != 2037L || any(invalid_key(sample$sample_id)) ||
    anyDuplicated(sample$sample_id)) {
  stop("Sample primary-key or row-count check failed: expected 2,037 unique UUIDs.")
}
if (any(invalid_key(sample$patient_id)) ||
    any(!(sample$patient_id %in% patient$patient_id))) {
  stop("Sample -> Patient foreign-key check failed.")
}
if (length(unique(sample$patient_id)) != 633L ||
    !setequal(sample$patient_id, patient$patient_id)) {
  stop("Sample cohort does not cover the 633 clinical cases as previously validated.")
}
if (!identical(sample$sample_id, sample$sample_sample_id)) {
  stop("Sample UUID differs from preserved source field.")
}

# Verify source project and case barcode against the parent Patient record.
parent_idx <- match(sample$patient_id, patient$patient_id)
project_mismatch <- (raw_sample[["project.project_id"]] !=
                     patient$project_project_id[parent_idx])
barcode_mismatch <- (raw_sample[["cases.submitter_id"]] !=
                     patient$case_submitter_id[parent_idx])
if (anyNA(project_mismatch) || any(project_mismatch) ||
    anyNA(barcode_mismatch) || any(barcode_mismatch)) {
  stop("Sample project / case barcode disagrees with the parent Patient record.")
}

project_counts <- table(raw_sample[["project.project_id"]])
if (!identical(unname(as.integer(project_counts[c("TCGA-COAD", "TCGA-READ")])),
               c(1490L, 547L))) {
  stop("Sample project distribution changed from source baseline.")
}

# STEP 4: Produce aggregate-only reports (no patient IDs or barcodes).
# TSV column presence is uniform, but values may be explicit GDC missing
# markers. Do not equate the marker "'--" with a negative clinical finding.
coverage <- do.call(rbind, lapply(source_sample_fields, function(field) {
  values <- raw_sample[[field]]
  data.frame(
    source_column = field,
    output_column = paste0("sample_", sub("^samples\\.", "", field)),
    records = length(values),
    gdc_missing_marker = sum(values == "'--", na.rm = TRUE),
    blank_values = sum(is.na(values) | values == ""),
    stringsAsFactors = FALSE
  )
}))

sample_types <- as.data.frame(
  table(sample$sample_sample_type, useNA = "ifany"),
  stringsAsFactors = FALSE
)
names(sample_types) <- c("sample_type", "records")

validation <- data.frame(
  metric = c(
    "sample_rows", "unique_sample_ids", "missing_sample_ids",
    "orphan_sample_patient_links", "patients_with_samples",
    "patients_without_samples", "sample_project_mismatches",
    "sample_case_barcode_mismatches", "sample_direct_source_fields",
    "coad_sample_rows", "read_sample_rows"
  ),
  value = c(
    nrow(sample), length(unique(sample$sample_id)),
    sum(invalid_key(sample$sample_id)),
    sum(!(sample$patient_id %in% patient$patient_id)),
    length(unique(sample$patient_id)),
    nrow(patient) - length(unique(sample$patient_id)),
    sum(project_mismatch), sum(barcode_mismatch),
    length(source_sample_fields),
    unname(project_counts["TCGA-COAD"]),
    unname(project_counts["TCGA-READ"])
  ),
  stringsAsFactors = FALSE
)

# STEP 5: Save only after all checks have passed.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sample, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)
utils::write.csv(sample_types, sample_types_file, row.names = FALSE)

cat("\nPhase 5 Sample extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Input source: biospecimen/sample.tsv (not the clinical JSON).\n")
cat("The GDC missing marker remains unchanged in source fields.\n")
cat("This is NOT yet an OMOP SPECIMEN table.\n")
