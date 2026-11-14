# TCGA COAD/READ -> OMOP | Phase 5: harmonised Analyte extraction
# Run from repository root:
#   source("R/03_harmonisation/12_build_analyte.R")
# Source: data/raw/gdc_r46/biospecimen/analyte.tsv (GDC Release 46.0)
# Grain: one row per analytes.analyte_id (GDC analyte UUID)
# Relationship: analyte.portion_id -> portion.portion_id
# No direct patient_id or sample_id in the harmonised Analyte table.
# All analytes.* fields remain source text; "'--" is NOT recoded here.
# Not an OMOP table.

source_file <- "data/raw/gdc_r46/biospecimen/analyte.tsv"
portion_file <- "data/processed/harmonised/portion.csv"
sample_file <- "data/processed/harmonised/sample.csv"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/analyte.csv"
validation_file <- "results/validation/12_analyte_validation.csv"
coverage_file <- "results/validation/12_analyte_field_coverage.csv"

for (f in c(source_file, portion_file, sample_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input: ", f, ". Run from project root.")
}

read_table <- function(f, sep) {
  utils::read.table(f, sep = sep, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character", na.strings = character(),
    quote = '"', comment.char = "", fileEncoding = "UTF-8-BOM")
}
raw_analyte <- read_table(source_file, "\t")
portion <- read_table(portion_file, ",")
sample <- read_table(sample_file, ",")
patient <- read_table(patient_file, ",")

require_fields <- function(x, fields, file_label) {
  missing <- setdiff(fields, names(x))
  if (length(missing)) stop(file_label, " missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(names(x))) stop(file_label, " has duplicate column names.")
}
require_fields(raw_analyte, c(
  "project.project_id", "cases.case_id", "cases.submitter_id",
  "samples.sample_id", "samples.submitter_id", "portions.portion_id",
  "portions.submitter_id", "analytes.analyte_id"), "analyte.tsv")
require_fields(portion, c("portion_id", "sample_id", "portion_submitter_id"), "portion.csv")
require_fields(sample, c("sample_id", "patient_id", "sample_submitter_id"), "sample.csv")
require_fields(patient, c("patient_id", "case_submitter_id", "project_project_id"), "patient.csv")

invalid_key <- function(x) is.na(x) | !nzchar(trimws(x)) | x == "'--"
check_parent <- function(x, key, expected_n, label) {
  if (nrow(x) != expected_n || any(invalid_key(x[[key]])) || anyDuplicated(x[[key]])) {
    stop(label, " does not match its validated primary-key baseline.")
  }
}
check_parent(portion, "portion_id", 3217L, "Portion")
check_parent(sample, "sample_id", 2037L, "Sample")
check_parent(patient, "patient_id", 633L, "Patient")

# STEP 1: Source-faithful extraction, with only direct parent FK.
analyte_fields <- grep("^analytes\\.", names(raw_analyte), value = TRUE)
if (nrow(raw_analyte) != 4436L || ncol(raw_analyte) != 23L ||
    length(analyte_fields) != 16L) {
  stop("Analyte source count/structure changed: got ", nrow(raw_analyte),
       " rows, ", ncol(raw_analyte), " columns, ", length(analyte_fields), " analytes.* fields")
}
analyte <- data.frame(
  analyte_id = raw_analyte[["analytes.analyte_id"]],
  portion_id = raw_analyte[["portions.portion_id"]],
  stringsAsFactors = FALSE, check.names = FALSE
)
for (field in analyte_fields) {
  new_name <- paste0("analyte_", sub("^analytes\\.", "", field))
  if (new_name %in% names(analyte)) stop("Duplicate output field: ", new_name)
  analyte[[new_name]] <- raw_analyte[[field]]
}

# STEP 2: Validate uniqueness and lineage against all verified parent tables.
if (any(invalid_key(analyte$analyte_id)) || anyDuplicated(analyte$analyte_id)) {
  stop("Analyte IDs must be nonmissing and unique.")
}
if (!identical(analyte$analyte_id, analyte$analyte_analyte_id)) {
  stop("Analyte UUID changed during extraction.")
}
if (any(invalid_key(analyte$portion_id))) stop("Invalid parent portion_id.")
portion_idx <- match(analyte$portion_id, portion$portion_id)
if (anyNA(portion_idx)) stop("Orphan Analyte -> Portion links.")
sample_idx <- match(portion$sample_id[portion_idx], sample$sample_id)
if (anyNA(sample_idx)) stop("Parent Portion -> Sample link broken.")
patient_idx <- match(sample$patient_id[sample_idx], patient$patient_id)
if (anyNA(patient_idx)) stop("Parent Sample -> Patient link broken.")

mismatch <- function(a, b) {
  if (anyNA(a) || anyNA(b)) stop("Unexpected NA in checked source linkage.")
  a != b
}
checks <- list(
  portion_barcode = mismatch(raw_analyte[["portions.submitter_id"]], portion$portion_submitter_id[portion_idx]),
  sample_id = mismatch(raw_analyte[["samples.sample_id"]], sample$sample_id[sample_idx]),
  sample_barcode = mismatch(raw_analyte[["samples.submitter_id"]], sample$sample_submitter_id[sample_idx]),
  case_id = mismatch(raw_analyte[["cases.case_id"]], patient$patient_id[patient_idx]),
  case_barcode = mismatch(raw_analyte[["cases.submitter_id"]], patient$case_submitter_id[patient_idx]),
  project = mismatch(raw_analyte[["project.project_id"]], patient$project_project_id[patient_idx])
)
if (any(vapply(checks, any, logical(1)))) {
  stop("Analyte ancestry mismatch: ", paste(names(checks)[vapply(checks, any, logical(1))], collapse = ", "))
}
projects <- table(raw_analyte[["project.project_id"]])
if (!identical(unname(as.integer(projects[c("TCGA-COAD", "TCGA-READ")])), c(3198L, 1238L))) {
  stop("Analyte COAD/READ counts differ from source baseline.")
}

# STEP 3: Aggregate validation and field-coverage outputs.
coverage <- do.call(rbind, lapply(analyte_fields, function(field) {
  x <- raw_analyte[[field]]
  data.frame(
    source_column = field,
    output_column = paste0("analyte_", sub("^analytes\\.", "", field)),
    records = length(x),
    gdc_missing_marker = sum(x == "'--", na.rm = TRUE),
    blank_values = sum(is.na(x) | x == ""),
    stringsAsFactors = FALSE
  )
}))
validation <- data.frame(
  metric = c("analyte_rows", "unique_analyte_ids", "missing_analyte_ids",
    "orphan_analyte_portion_links", "orphan_analyte_sample_links", "orphan_analyte_patient_links",
    "portion_barcode_mismatches", "sample_id_mismatches", "sample_barcode_mismatches",
    "case_id_mismatches", "case_barcode_mismatches", "project_mismatches",
    "portions_with_analytes", "samples_with_analytes", "patients_with_analytes",
    "analyte_direct_source_fields", "coad_analyte_rows", "read_analyte_rows"),
  value = c(nrow(analyte), length(unique(analyte$analyte_id)), sum(invalid_key(analyte$analyte_id)),
    sum(is.na(portion_idx)), sum(is.na(sample_idx)), sum(is.na(patient_idx)),
    sum(checks$portion_barcode), sum(checks$sample_id), sum(checks$sample_barcode),
    sum(checks$case_id), sum(checks$case_barcode), sum(checks$project),
    length(unique(analyte$portion_id)), length(unique(portion$sample_id[portion_idx])),
    length(unique(patient$patient_id[patient_idx])), length(analyte_fields),
    unname(projects["TCGA-COAD"]), unname(projects["TCGA-READ"])),
  stringsAsFactors = FALSE
)

# STEP 4: Save only if checks have all passed.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(analyte, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)
cat("\nPhase 5 Analyte extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Source: biospecimen/analyte.tsv; Analyte links to Portion.\n")
cat("Some Portions lack analyte records; absence is not a negative laboratory finding.\n")
cat("Unconverted source values and GDC missing markers remain intact.\n")
cat("This is NOT yet an OMOP SPECIMEN table.\n")
