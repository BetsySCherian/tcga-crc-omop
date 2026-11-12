# TCGA COAD/READ -> OMOP | Phase 5: harmonised Slide and Slide--Portion
# Run from repository root:
#   source("R/03_harmonisation/14_build_slide_and_slide_portion.R")
#
# Source: GDC Release 46.0, biospecimen/slide.tsv.
# Resolved Phase 4 decision:
#   - slide: exactly one row per GDC slide UUID (1998 rows)
#   - slide_portion: exactly one row per (slide_id, portion_id) (2035 rows)
#   - 35 slide IDs relate to multiple portions; DO NOT remove these relationships
#   - each repeated slide ID has consistent slide attributes and belongs to one sample
# Foreign keys: slide_portion.slide_id -> slide.slide_id;
#               slide_portion.portion_id -> portion.portion_id.
# Parent chain for sample and patient comes via portion -> sample -> patient.
# Every slides.* field is retained as verbatim character text in slide.
# Uninterpreted GDC missing marker "'--" is preserved.
# This is NOT an OMOP SPECIMEN table.

source_file <- "data/raw/gdc_r46/biospecimen/slide.tsv"
portion_file <- "data/processed/harmonised/portion.csv"
sample_file <- "data/processed/harmonised/sample.csv"
patient_file <- "data/processed/harmonised/patient.csv"
slide_file <- "data/processed/harmonised/slide.csv"
slide_portion_file <- "data/processed/harmonised/slide_portion.csv"
validation_file <- "results/validation/14_slide_portion_validation.csv"
coverage_file <- "results/validation/14_slide_field_coverage.csv"
frequency_file <- "results/validation/14_slide_portion_frequency.csv"

for (file in c(source_file, portion_file, sample_file, patient_file)) {
  if (!file.exists(file)) stop("Missing required input: ", file, "\nRun from repository root.")
}

# STEP 1: Read all identifiers and original source fields without conversion.
read_table <- function(file, separator) {
  utils::read.table(
    file, header = TRUE, sep = separator, quote = '"', comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = character(), fileEncoding = "UTF-8-BOM"
  )
}
raw_slide <- read_table(source_file, "\t")
portion <- read_table(portion_file, ",")
sample <- read_table(sample_file, ",")
patient <- read_table(patient_file, ",")

require_fields <- function(tbl, fields, label) {
  missing <- setdiff(fields, names(tbl))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(names(tbl))) stop(label, " contains duplicate column names.")
}
require_fields(raw_slide, c(
  "project.project_id", "cases.case_id", "cases.submitter_id",
  "samples.sample_id", "samples.submitter_id",
  "portions.portion_id", "portions.submitter_id",
  "slides.slide_id", "slides.submitter_id"
), "slide.tsv")
require_fields(portion, c("portion_id", "sample_id", "portion__submitter_id"), "portion.csv")
require_fields(sample, c("sample_id", "patient_id", "sample__submitter_id"), "sample.csv")
require_fields(patient, c("patient_id", "case__submitter_id", "project__project_id"), "patient.csv")

invalid_key <- function(x) is.na(x) | !nzchar(trimws(x)) | x == "'--"
check_parent <- function(tbl, col, expected, label) {
  if (nrow(tbl) != expected || any(invalid_key(tbl[[col]])) || anyDuplicated(tbl[[col]])) {
    stop(label, " parent table fails row-count or key checks.")
  }
}
check_parent(portion, "portion_id", 3217L, "Portion")
check_parent(sample, "sample_id", 2037L, "Sample")
check_parent(patient, "patient_id", 633L, "Patient")

slide_fields <- grep("^slides\\.", names(raw_slide), value = TRUE)
if (nrow(raw_slide) != 2035L || ncol(raw_slide) != 31L || length(slide_fields) != 24L) {
  stop("Source slide.tsv does not match verified GDC R46 baseline: ", nrow(raw_slide),
       " rows, ", ncol(raw_slide), " columns, ", length(slide_fields), " slide fields.")
}

# STEP 2: Validate original slide--portion associations and full source ancestry.
source_slide_ids <- raw_slide[["slides.slide_id"]]
source_portion_ids <- raw_slide[["portions.portion_id"]]
if (any(invalid_key(source_slide_ids)) || any(invalid_key(source_portion_ids))) {
  stop("Missing slide ID or portion ID in slide.tsv.")
}

pair_key <- paste(source_slide_ids, source_portion_ids, sep = "::")
if (anyDuplicated(pair_key)) stop("Duplicate slide--portion pairs; cannot use composite primary key.")

portion_idx <- match(source_portion_ids, portion$portion_id)
if (anyNA(portion_idx)) stop("Slide--Portion -> Portion orphan links.")
sample_idx <- match(portion$sample_id[portion_idx], sample$sample_id)
if (anyNA(sample_idx)) stop("Slide--Portion -> Sample orphan links.")
patient_idx <- match(sample$patient_id[sample_idx], patient$patient_id)
if (anyNA(patient_idx)) stop("Slide--Portion -> Patient orphan links.")

compare <- function(left, right, label) {
  if (length(left) != length(right) || anyNA(left) || anyNA(right)) {
    stop("Unexpected missing value in relationship check: ", label)
  }
  left != right
}
checks <- list(
  portion_barcode = compare(raw_slide[["portions.submitter_id"]],
    portion$portion__submitter_id[portion_idx], "portion barcode"),
  sample_id = compare(raw_slide[["samples.sample_id"]],
    sample$sample_id[sample_idx], "sample UUID"),
  sample_barcode = compare(raw_slide[["samples.submitter_id"]],
    sample$sample__submitter_id[sample_idx], "sample barcode"),
  case_id = compare(raw_slide[["cases.case_id"]],
    patient$patient_id[patient_idx], "case UUID"),
  case_barcode = compare(raw_slide[["cases.submitter_id"]],
    patient$case__submitter_id[patient_idx], "case barcode"),
  project = compare(raw_slide[["project.project_id"]],
    patient$project__project_id[patient_idx], "project")
)
failed_checks <- names(checks)[vapply(checks, any, logical(1))]
if (length(failed_checks)) stop("Inconsistent slide ancestry: ", paste(failed_checks, collapse = ", "))

# STEP 3: Confirm duplicate slide UUIDs have identical slide-level attributes.
# Selecting the first source row per slide is lossless only if all its slide
# attributes are identical across the other rows for that same slide UUID.
first_indices <- !duplicated(source_slide_ids)
first_index_for_row <- match(source_slide_ids, source_slide_ids[first_indices])
first_source_rows <- which(first_indices)
slide_attribute_mismatch <- rep(FALSE, nrow(raw_slide))
for (field in slide_fields) {
  values <- raw_slide[[field]]
  slide_attribute_mismatch <- slide_attribute_mismatch |
    (values != values[first_source_rows][first_index_for_row])
}
if (anyNA(slide_attribute_mismatch) || any(slide_attribute_mismatch)) {
  stop("One slide UUID has conflicting slides.* source attributes across portions.")
}

# Confirm one sample per distinct slide (but potentially multiple portions).
slide_sample_conflict <- raw_slide[["samples.sample_id"]] !=
  raw_slide[["samples.sample_id"]][first_source_rows][first_index_for_row]
if (anyNA(slide_sample_conflict) || any(slide_sample_conflict)) {
  stop("A slide UUID appears under more than one sample; review source model.")
}

# STEP 4: Materialise both entities, with separate grains.
# One record per distinct physical slide UUID; all 24 slide attributes retained.
slide <- data.frame(slide_id = source_slide_ids[first_indices],
                    stringsAsFactors = FALSE, check.names = FALSE)
for (field in slide_fields) {
  output_name <- paste0("slide__", sub("^slides\\.", "", field))
  if (output_name %in% names(slide)) stop("Slide output column collision: ", output_name)
  slide[[output_name]] <- raw_slide[[field]][first_indices]
}
# One record per observed slide / portion relationship.
slide_portion <- data.frame(
  slide_id = source_slide_ids,
  portion_id = source_portion_ids,
  source_tsv_row = seq_len(nrow(raw_slide)),
  stringsAsFactors = FALSE, check.names = FALSE
)

if (nrow(slide) != 1998L || any(invalid_key(slide$slide_id)) ||
    anyDuplicated(slide$slide_id)) stop("Slide table primary-key/count checks failed.")
if (!identical(slide$slide_id, slide$slide__slide_id)) {
  stop("Slide ID does not match preserved slides.slide_id field.")
}
if (nrow(slide_portion) != 2035L || anyDuplicated(slide_portion[c("slide_id", "portion_id")])) {
  stop("Slide--Portion table count/composite-key check failed.")
}
if (anyNA(match(slide_portion$slide_id, slide$slide_id))) {
  stop("Slide--Portion -> Slide orphan links.")
}
if (anyNA(match(slide_portion$portion_id, portion$portion_id))) {
  stop("Slide--Portion -> Portion orphan links.")
}
if (!identical(slide_portion$slide_id, source_slide_ids) ||
    !identical(slide_portion$portion_id, source_portion_ids)) {
  stop("Source slide--portion associations were altered.")
}

frequency <- as.integer(table(source_slide_ids))
slide_frequency <- data.frame(
  portions_per_slide = as.integer(names(table(frequency))),
  number_of_slides = as.integer(table(frequency)),
  stringsAsFactors = FALSE
)
repeated_ids <- sum(frequency > 1L)
project_counts <- table(raw_slide[["project.project_id"]])
if (!identical(unname(as.integer(project_counts[c("TCGA-COAD", "TCGA-READ")])),
               c(1493L, 542L))) stop("Slide--Portion source project counts changed.")
if (repeated_ids != 35L ||
    !identical(unname(as.integer(table(frequency)[c("1", "2", "3")])),
               c(1963L, 33L, 2L))) stop("Repeated-slide pattern differs from resolved exploration.")
if (length(unique(raw_slide[["samples.sample_id"]])) != 1372L ||
    length(unique(raw_slide[["cases.case_id"]])) != 631L) {
  stop("Slide-associated sample / patient counts differ from source baseline.")
}

# STEP 5: Aggregate-only reports, with no UUIDs or patient barcodes.
coverage <- do.call(rbind, lapply(slide_fields, function(field) {
  values <- raw_slide[[field]][first_indices]
  data.frame(
    source_column = field,
    output_column = paste0("slide__", sub("^slides\\.", "", field)),
    slides = length(values),
    gdc_missing_marker = sum(values == "'--", na.rm = TRUE),
    blank_values = sum(is.na(values) | values == ""),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  metric = c(
    "slide_source_rows", "slide_rows", "unique_slide_ids", "missing_slide_ids",
    "slide_portion_rows", "unique_slide_portion_pairs",
    "repeated_slide_ids", "slide_attribute_conflicts", "slide_sample_conflicts",
    "orphan_slide_portion_slide_links", "orphan_slide_portion_portion_links",
    "orphan_slide_portion_sample_links", "orphan_slide_portion_patient_links",
    paste0(names(checks), "_mismatches"),
    "portions_with_slides", "samples_with_slides", "patients_with_slides",
    "slide_direct_source_fields", "coad_slide_portion_rows", "read_slide_portion_rows"
  ),
  value = c(
    nrow(raw_slide), nrow(slide), length(unique(slide$slide_id)),
    sum(invalid_key(slide$slide_id)), nrow(slide_portion),
    nrow(unique(slide_portion[c("slide_id", "portion_id")])),
    repeated_ids, sum(slide_attribute_mismatch), sum(slide_sample_conflict),
    sum(is.na(match(slide_portion$slide_id, slide$slide_id))),
    sum(is.na(portion_idx)), sum(is.na(sample_idx)), sum(is.na(patient_idx)),
    vapply(checks, function(v) as.integer(sum(v)), integer(1)),
    length(unique(source_portion_ids)), length(unique(raw_slide[["samples.sample_id"]])),
    length(unique(raw_slide[["cases.case_id"]])), length(slide_fields),
    unname(project_counts["TCGA-COAD"]), unname(project_counts["TCGA-READ"])
  ), stringsAsFactors = FALSE
)

# STEP 6: Write outputs only after all validations pass.
# Data outputs contain row-level source UUIDs; keep Git-ignored.
dir.create(dirname(slide_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(slide, slide_file, row.names = FALSE, na = "")
utils::write.csv(slide_portion, slide_portion_file, row.names = FALSE, na = "")
utils::write.csv(validation, validation_file, row.names = FALSE)
utils::write.csv(coverage, coverage_file, row.names = FALSE)
utils::write.csv(slide_frequency, frequency_file, row.names = FALSE)

cat("\nPhase 5 Slide and Slide--Portion extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nOutputs:\n  ", slide_file, "\n  ", slide_portion_file, "\n", sep = "")
cat("The 35 slides with multiple portion links are preserved without duplicate Slide entities.\n")
cat("These are NOT yet OMOP SPECIMEN tables.\n")
