# TCGA COAD/READ -> OMOP | Phase 5: harmonised Portion extraction
# Run from the repository root:
#   source("R/03_harmonisation/11_build_portion.R")
#
# Source: GDC Release 46.0 biospecimen/portion.tsv.
# Grain: one row per GDC portions.portion_id (portion UUID).
# Relationship: portion.sample_id -> sample.sample_id.
# Do not repeat patient_id on Portion: derive it via Sample when needed.
# All portions.* values are preserved as text, including the GDC missing
# marker "'--" and the raw creation_datetime (Unix timestamp where present).
# This is NOT an OMOP SPECIMEN table.

source_file <- "data/raw/gdc_r46/biospecimen/portion.tsv"
sample_file <- "data/processed/harmonised/sample.csv"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/portion.csv"
validation_file <- "results/validation/11_portion_validation.csv"
coverage_file <- "results/validation/11_portion_field_coverage.csv"

for (f in c(source_file, sample_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input: ", f, "\nRun from the project root.")
}

# STEP 1. Import source and previously validated parent tables as text.
# Do not interpret GDC missing markers as negative results or convert dates.
raw_portion <- utils::read.delim(
  source_file, sep = "\t", header = TRUE,
  check.names = FALSE, stringsAsFactors = FALSE,
  colClasses = "character", na.strings = character(),
  quote = "\"", comment.char = "", fileEncoding = "UTF-8-BOM"
)
sample <- utils::read.csv(
  sample_file, check.names = FALSE, stringsAsFactors = FALSE,
  colClasses = "character", na.strings = character()
)
patient <- utils::read.csv(
  patient_file, check.names = FALSE, stringsAsFactors = FALSE,
  colClasses = "character", na.strings = character()
)

source_required <- c(
  "project.project_id", "cases.case_id", "cases.submitter_id",
  "samples.sample_id", "samples.submitter_id", "portions.portion_id"
)
parent_sample_required <- c("sample_id", "patient_id", "sample_submitter_id")
parent_patient_required <- c("patient_id", "case_submitter_id", "project_project_id")
for (spec in list(
  list(x = raw_portion, fields = source_required, label = "portion.tsv"),
  list(x = sample, fields = parent_sample_required, label = "sample.csv"),
  list(x = patient, fields = parent_patient_required, label = "patient.csv")
)) {
  missing <- setdiff(spec$fields, names(spec$x))
  if (length(missing)) stop(spec$label, " is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(names(spec$x))) stop("Duplicate column names in ", spec$label)
}

invalid_key <- function(x) is.na(x) | !nzchar(trimws(x)) | x == "'--"
if (nrow(sample) != 2037L || any(invalid_key(sample$sample_id)) ||
    anyDuplicated(sample$sample_id)) stop("Sample parent table fails its baseline primary-key checks.")
if (nrow(patient) != 633L || any(invalid_key(patient$patient_id)) ||
    anyDuplicated(patient$patient_id)) stop("Patient parent table fails its baseline primary-key checks.")

# STEP 2. Create one row per portion and keep every portions.* source attribute.
portion_fields <- grep("^portions\\.", names(raw_portion), value = TRUE)
if (nrow(raw_portion) != 3217L || length(portion_fields) != 7L ||
    ncol(raw_portion) != 12L) {
  stop("Portion source structure/count changed: expected 3217 rows, 7 portions.* fields, 12 TSV columns; observed ",
       nrow(raw_portion), ", ", length(portion_fields), ", ", ncol(raw_portion))
}

portion <- data.frame(
  portion_id = raw_portion[["portions.portion_id"]],
  sample_id = raw_portion[["samples.sample_id"]],
  stringsAsFactors = FALSE, check.names = FALSE
)
for (field in portion_fields) {
  output_name <- paste0("portion_", sub("^portions\\.", "", field))
  if (output_name %in% names(portion)) stop("Output column name collision: ", output_name)
  portion[[output_name]] <- raw_portion[[field]]
}

# STEP 3. Validate primary keys, foreign keys, and hierarchy coherence.
if (any(invalid_key(portion$portion_id)) || anyDuplicated(portion$portion_id)) {
  stop("Portion primary key must be nonmissing and unique for 3217 records.")
}
if (!identical(portion$portion_id, portion$portion_portion_id)) {
  stop("Preserved portions.portion_id differs from portion_id.")
}
if (any(invalid_key(portion$sample_id))) stop("Missing/invalid sample link on Portion.")

sample_idx <- match(portion$sample_id, sample$sample_id)
if (anyNA(sample_idx)) stop("Orphan Portion -> Sample links detected.")
case_mismatch <- raw_portion[["cases.case_id"]] != sample$patient_id[sample_idx]
sample_barcode_mismatch <- raw_portion[["samples.submitter_id"]] != sample$sample_submitter_id[sample_idx]
if (anyNA(case_mismatch) || any(case_mismatch) ||
    anyNA(sample_barcode_mismatch) || any(sample_barcode_mismatch)) {
  stop("Portion source case/sample barcode conflicts with parent Sample.")
}

patient_idx <- match(sample$patient_id[sample_idx], patient$patient_id)
if (anyNA(patient_idx)) stop("Portion's parent Sample points to a missing Patient.")
case_barcode_mismatch <- raw_portion[["cases.submitter_id"]] !=
  patient$case_submitter_id[patient_idx]
project_mismatch <- raw_portion[["project.project_id"]] !=
  patient$project_project_id[patient_idx]
if (anyNA(case_barcode_mismatch) || any(case_barcode_mismatch) ||
    anyNA(project_mismatch) || any(project_mismatch)) {
  stop("Portion case barcode/project conflicts with parent Patient.")
}

if (!setequal(unique(portion$sample_id), sample$sample_id)) {
  stop("Portion coverage differs from expected: all 2037 Samples should have portions.")
}
project_counts <- table(raw_portion[["project.project_id"]])
if (!identical(unname(as.integer(project_counts[c("TCGA-COAD", "TCGA-READ")])),
               c(2357L, 860L))) stop("Portion project counts differ from verified baseline.")

# STEP 4. Aggregate-only validation and source-field coverage reports.
coverage <- do.call(rbind, lapply(portion_fields, function(field) {
  values <- raw_portion[[field]]
  data.frame(
    source_column = field,
    output_column = paste0("portion_", sub("^portions\\.", "", field)),
    records = length(values),
    gdc_missing_marker = sum(values == "'--", na.rm = TRUE),
    blank_values = sum(is.na(values) | values == ""),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "portion_rows", "unique_portion_ids", "missing_portion_ids",
    "orphan_portion_sample_links", "orphan_portion_patient_links",
    "portion_case_id_mismatches", "portion_sample_barcode_mismatches",
    "portion_case_barcode_mismatches", "portion_project_mismatches",
    "samples_with_portions", "patients_with_portions",
    "portion_direct_source_fields", "coad_portion_rows", "read_portion_rows"
  ),
  value = c(
    nrow(portion), length(unique(portion$portion_id)), sum(invalid_key(portion$portion_id)),
    sum(is.na(sample_idx)), sum(is.na(patient_idx)),
    sum(case_mismatch), sum(sample_barcode_mismatch),
    sum(case_barcode_mismatch), sum(project_mismatch),
    length(unique(portion$sample_id)), length(unique(sample$patient_id[sample_idx])),
    length(portion_fields), unname(project_counts["TCGA-COAD"]),
    unname(project_counts["TCGA-READ"])
  ),
  stringsAsFactors = FALSE
)

# STEP 5. Write outputs only after validation succeeds.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(portion, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 Portion extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Input source: biospecimen/portion.tsv, linked to the existing Sample table.\n")
cat("Creation timestamps and GDC missing markers remain unconverted source values.\n")
cat("This is NOT yet an OMOP SPECIMEN table.\n")
