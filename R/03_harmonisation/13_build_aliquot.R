# TCGA COAD/READ -> OMOP | Phase 5: harmonised Aliquot extraction
# Run from the repository root:
#   source("R/03_harmonisation/13_build_aliquot.R")
#
# Source: GDC Release 46.0 biospecimen/aliquot.tsv.
# Grain: one row per GDC aliquots.aliquot_id (unique aliquot UUID).
# Relationship: aliquot.analyte_id -> analyte.analyte_id.
# Do not duplicate patient_id/sample_id/portion_id in the Aliquot entity:
# those relationships are derived through Analyte -> Portion -> Sample -> Patient.
# Preserve all aliquots.* fields verbatim as character text, including the
# GDC missing-value marker "'--". This is NOT an OMOP SPECIMEN table.

source_file <- "data/raw/gdc_r46/biospecimen/aliquot.tsv"
analyte_file <- "data/processed/harmonised/analyte.csv"
portion_file <- "data/processed/harmonised/portion.csv"
sample_file <- "data/processed/harmonised/sample.csv"
patient_file <- "data/processed/harmonised/patient.csv"
out_file <- "data/processed/harmonised/aliquot.csv"
validation_file <- "results/validation/13_aliquot_validation.csv"
coverage_file <- "results/validation/13_aliquot_field_coverage.csv"

for (f in c(source_file, analyte_file, portion_file, sample_file, patient_file)) {
  if (!file.exists(f)) stop("Missing input: ", f, "\nRun from the project root.")
}

# STEP 1. Import as text to retain identifiers and exact GDC source values.
read_table <- function(file, sep) {
  utils::read.table(
    file, header = TRUE, sep = sep, quote = '"', comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character(),
    fileEncoding = "UTF-8-BOM"
  )
}
raw_aliquot <- read_table(source_file, "\t")
analyte <- read_table(analyte_file, ",")
portion <- read_table(portion_file, ",")
sample <- read_table(sample_file, ",")
patient <- read_table(patient_file, ",")

require_fields <- function(data, fields, label) {
  missing <- setdiff(fields, names(data))
  if (length(missing)) stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(names(data))) stop(label, " has duplicate column names.")
}
require_fields(raw_aliquot, c(
  "project.project_id", "cases.case_id", "cases.submitter_id",
  "samples.sample_id", "samples.submitter_id",
  "portions.portion_id", "portions.submitter_id",
  "analytes.analyte_id", "analytes.submitter_id",
  "aliquots.aliquot_id"
), "aliquot.tsv")
require_fields(analyte, c("analyte_id", "portion_id", "analyte__submitter_id"), "analyte.csv")
require_fields(portion, c("portion_id", "sample_id", "portion__submitter_id"), "portion.csv")
require_fields(sample, c("sample_id", "patient_id", "sample__submitter_id"), "sample.csv")
require_fields(patient, c("patient_id", "case__submitter_id", "project__project_id"), "patient.csv")

invalid_key <- function(x) is.na(x) | !nzchar(trimws(x)) | x == "'--"
check_parent <- function(data, key, expected_rows, label) {
  if (nrow(data) != expected_rows || any(invalid_key(data[[key]])) ||
      anyDuplicated(data[[key]])) {
    stop(label, " does not match the validated parent primary-key baseline.")
  }
}
check_parent(analyte, "analyte_id", 4436L, "Analyte")
check_parent(portion, "portion_id", 3217L, "Portion")
check_parent(sample, "sample_id", 2037L, "Sample")
check_parent(patient, "patient_id", 633L, "Patient")

# STEP 2. Construct the source-faithful Aliquot entity.
aliquot_fields <- grep("^aliquots\\.", names(raw_aliquot), value = TRUE)
if (nrow(raw_aliquot) != 12574L || ncol(raw_aliquot) != 26L ||
    length(aliquot_fields) != 17L) {
  stop("Aliquot source differs from verified GDC R46 baseline: ",
       nrow(raw_aliquot), " rows, ", ncol(raw_aliquot), " columns, ",
       length(aliquot_fields), " aliquots.* fields.")
}
aliquot <- data.frame(
  aliquot_id = raw_aliquot[["aliquots.aliquot_id"]],
  analyte_id = raw_aliquot[["analytes.analyte_id"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (field in aliquot_fields) {
  new_name <- paste0("aliquot__", sub("^aliquots\\.", "", field))
  if (new_name %in% names(aliquot)) stop("Aliquot output field collision: ", new_name)
  aliquot[[new_name]] <- raw_aliquot[[field]]
}

# STEP 3. Validate row grain and full chain of parent relationships.
if (any(invalid_key(aliquot$aliquot_id)) || anyDuplicated(aliquot$aliquot_id)) {
  stop("Aliquot primary keys are not all nonmissing and unique.")
}
if (!identical(aliquot$aliquot_id, aliquot$aliquot__aliquot_id)) {
  stop("Aliquot ID differs from its original source field.")
}
if (any(invalid_key(aliquot$analyte_id))) stop("Invalid Analyte foreign key in Aliquot.")

analyte_idx <- match(aliquot$analyte_id, analyte$analyte_id)
if (anyNA(analyte_idx)) stop("Orphan Aliquot -> Analyte links.")
portion_idx <- match(analyte$portion_id[analyte_idx], portion$portion_id)
if (anyNA(portion_idx)) stop("Orphan Analyte -> Portion ancestor links.")
sample_idx <- match(portion$sample_id[portion_idx], sample$sample_id)
if (anyNA(sample_idx)) stop("Orphan Portion -> Sample ancestor links.")
patient_idx <- match(sample$patient_id[sample_idx], patient$patient_id)
if (anyNA(patient_idx)) stop("Orphan Sample -> Patient ancestor links.")

# Cross-check every denormalised ancestor identifier/barcode in source TSV.
compare_values <- function(original, from_parent, label) {
  if (length(original) != length(from_parent) ||
      anyNA(original) || anyNA(from_parent)) {
    stop("Unexpected missing value in checked lineage: ", label)
  }
  original != from_parent
}
checks <- list(
  analyte_barcode = compare_values(raw_aliquot[["analytes.submitter_id"]],
                                   analyte$analyte__submitter_id[analyte_idx], "analyte barcode"),
  portion_id = compare_values(raw_aliquot[["portions.portion_id"]],
                              portion$portion_id[portion_idx], "portion UUID"),
  portion_barcode = compare_values(raw_aliquot[["portions.submitter_id"]],
                                   portion$portion__submitter_id[portion_idx], "portion barcode"),
  sample_id = compare_values(raw_aliquot[["samples.sample_id"]],
                             sample$sample_id[sample_idx], "sample UUID"),
  sample_barcode = compare_values(raw_aliquot[["samples.submitter_id"]],
                                  sample$sample__submitter_id[sample_idx], "sample barcode"),
  case_id = compare_values(raw_aliquot[["cases.case_id"]],
                           patient$patient_id[patient_idx], "case UUID"),
  case_barcode = compare_values(raw_aliquot[["cases.submitter_id"]],
                                patient$case__submitter_id[patient_idx], "case barcode"),
  project = compare_values(raw_aliquot[["project.project_id"]],
                           patient$project__project_id[patient_idx], "project")
)
failures <- names(checks)[vapply(checks, any, logical(1))]
if (length(failures)) stop("Aliquot source/parent mismatches in: ", paste(failures, collapse = ", "))

project_counts <- table(raw_aliquot[["project.project_id"]])
if (!identical(unname(as.integer(project_counts[c("TCGA-COAD", "TCGA-READ")])),
               c(9159L, 3415L))) {
  stop("Aliquot project counts differ from the GDC R46 baseline.")
}
if (!setequal(unique(aliquot$analyte_id), analyte$analyte_id)) {
  stop("Expected every one of the 4436 analytes to have an aliquot.")
}

# STEP 4. Create aggregate-only reports, with no patient IDs or barcodes.
coverage <- do.call(rbind, lapply(aliquot_fields, function(field) {
  v <- raw_aliquot[[field]]
  data.frame(
    source_column = field,
    output_column = paste0("aliquot__", sub("^aliquots\\.", "", field)),
    records = length(v),
    gdc_missing_marker = sum(v == "'--", na.rm = TRUE),
    blank_values = sum(is.na(v) | v == ""),
    stringsAsFactors = FALSE
  )
}))
validation <- data.frame(
  metric = c(
    "aliquot_rows", "unique_aliquot_ids", "missing_aliquot_ids",
    "orphan_aliquot_analyte_links", "orphan_aliquot_portion_ancestry",
    "orphan_aliquot_sample_ancestry", "orphan_aliquot_patient_ancestry",
    paste0(names(checks), "_mismatches"),
    "analytes_with_aliquots", "portions_with_aliquots", "samples_with_aliquots",
    "patients_with_aliquots", "aliquot_direct_source_fields",
    "coad_aliquot_rows", "read_aliquot_rows"
  ),
  value = c(
    nrow(aliquot), length(unique(aliquot$aliquot_id)), sum(invalid_key(aliquot$aliquot_id)),
    sum(is.na(analyte_idx)), sum(is.na(portion_idx)),
    sum(is.na(sample_idx)), sum(is.na(patient_idx)),
    vapply(checks, function(x) as.integer(sum(x)), integer(1)),
    length(unique(aliquot$analyte_id)), length(unique(analyte$portion_id[analyte_idx])),
    length(unique(portion$sample_id[portion_idx])), length(unique(patient$patient_id[patient_idx])),
    length(aliquot_fields), unname(project_counts["TCGA-COAD"]),
    unname(project_counts["TCGA-READ"])
  ),
  stringsAsFactors = FALSE
)

# STEP 5. Write outputs only if all checks pass.
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(aliquot, out_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)

cat("\nPhase 5 Aliquot extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutput: ", out_file, "\n", sep = "")
cat("Source: biospecimen/aliquot.tsv. Aliquot links to Analyte.\n")
cat("All aliquots.* fields were preserved as source text.\n")
cat("This is NOT yet an OMOP SPECIMEN table.\n")
