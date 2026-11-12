# TCGA COAD/READ -> OMOP CDM v5.4 | Phase 5
# Source-preserving GDC clinical extraction + patient table + schema validation
# Run at repository root: source("R/03_harmonisation/01_build_patient.R")
# Input: GDC Release 46.0 clinical JSON + existing exploration schema CSV
# This does NOT map to OMOP PERSON; no synthetic dates or concepts are generated.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install the required package first: install.packages('jsonlite')")
}

source_file <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
schema_file <- "results/exploration/01_gdc_json_schema.csv"
out_dir <- "data/processed/harmonised"
report_dir <- "results/validation"

if (!file.exists(source_file)) stop("Missing clinical JSON: ", source_file,
                                    "\nSet your RStudio working directory to the repository root.")
if (!file.exists(schema_file)) stop("Missing existing schema inventory: ", schema_file,
                                    "\nCopy exploration/01_gdc_json_schema.csv into results/exploration/.")

cases <- jsonlite::fromJSON(source_file, simplifyVector = FALSE)
if (!is.list(cases) || length(cases) != 633L) {
  stop("Expected 633 GDC cases, found ", length(cases), "; inspect source before continuing.")
}

# Extract all scalar keys at the case, demographic and project levels.
# Other nested data are not flattened into patient: they are retained below
# and will become separate linked tables in later scripts.
scalar_value <- function(obj, key) {
  if (!is.list(obj) || !(key %in% names(obj))) return(NA_character_)
  v <- obj[[key]]
  if (is.null(v) || is.list(v) || length(v) != 1L || is.na(v)) return(NA_character_)
  if (is.logical(v)) return(tolower(as.character(v)))
  as.character(v)
}
scalar_keys <- function(objs) {
  sort(unique(unlist(lapply(objs, function(obj) {
    if (!is.list(obj)) return(character())
    names(obj)[vapply(obj, function(v) !is.list(v) && length(v) == 1L,
                       logical(1))]
  }), use.names = FALSE)))
}
case_keys <- scalar_keys(cases)
demo_keys <- scalar_keys(lapply(cases, function(x) x$demographic))
project_keys <- scalar_keys(lapply(cases, function(x) x$project))

patient_rows <- lapply(cases, function(x) {
  c(patient_id = scalar_value(x, "case_id"),
    setNames(vapply(case_keys, function(k) scalar_value(x, k), character(1)),
             paste0("case__", case_keys)),
    setNames(vapply(demo_keys, function(k) scalar_value(x$demographic, k), character(1)),
             paste0("demographic__", demo_keys)),
    setNames(vapply(project_keys, function(k) scalar_value(x$project, k), character(1)),
             paste0("project__", project_keys)))
})
patient <- as.data.frame(do.call(rbind, patient_rows),
                         stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(patient) != 633L || anyNA(patient$patient_id) ||
    anyDuplicated(patient$patient_id)) {
  stop("Patient record count or primary key validation failed; no output saved.")
}
if (!"project__project_id" %in% names(patient)) stop("Project ID missing")
project_counts <- table(patient$project__project_id)
if (!identical(unname(as.integer(project_counts[c("TCGA-COAD", "TCGA-READ")])),
               c(461L, 172L))) {
  stop("Project counts differ from the verified 461 COAD / 172 READ baseline.")
}

# Retain every nested node with a reversible JSON Pointer (RFC 6901).
# For complete, byte-for-byte fidelity, copy and checksum the original JSON.
# The complete JSON snapshot is authoritative for distinguishing empty [] / {}.
escape_pointer <- function(s) gsub("/", "~1", gsub("~", "~0", s, fixed = TRUE), fixed = TRUE)
decode_pointer <- function(s) gsub("~0", "~", gsub("~1", "/", s, fixed = TRUE), fixed = TRUE)
walk_case <- function(root, id) {
  rows <- list()
  add <- function(path, typ, value, n) {
    rows[[length(rows) + 1L]] <<- data.frame(
      patient_id = id, json_pointer = path, node_type = typ,
      value_json = value, child_count = as.integer(n),
      stringsAsFactors = FALSE)
  }
  walk <- function(x, path) {
    if (is.null(x)) {
      add(path, "null", "null", 0L)
    } else if (is.list(x)) {
      obj <- !is.null(names(x))
      add(path, if (obj) "object" else "array", NA_character_, length(x))
      for (j in seq_along(x)) {
        token <- if (obj) escape_pointer(names(x)[j]) else as.character(j - 1L)
        walk(x[[j]], paste0(path, "/", token))
      }
    } else if (length(x) == 1L) {
      typ <- if (is.character(x)) "string" else if (is.logical(x)) "boolean" else "number"
      if (is.na(x)) add(path, "null", "null", 0L)
      else add(path, typ, as.character(jsonlite::toJSON(x, auto_unbox = TRUE,
                                                        digits = NA)), 0L)
    } else {
      stop("Unexpected atomic vector at JSON pointer: ", path)
    }
  }
  walk(root, "")
  do.call(rbind, rows)
}

attributes <- do.call(rbind, lapply(seq_along(cases), function(i) {
  walk_case(cases[[i]], patient$patient_id[i])
}))
rownames(attributes) <- NULL
if (anyDuplicated(paste(attributes$patient_id, attributes$json_pointer, sep = "\t"))) {
  stop("Duplicate patient / JSON-pointer combinations; no output saved.")
}

# Convert each observed JSON Pointer into the notation used by the
# PREVIOUSLY PRODUCED 01_gdc_json_schema.csv, e.g.:
# /diagnoses/0/treatments/1 -> root[].diagnoses[].treatments[]
pointer_to_schema <- function(ptr) {
  if (ptr == "") return("root[]")
  tokens <- strsplit(substring(ptr, 2L), "/", fixed = TRUE)[[1]]
  output <- "root[]"
  for (token in tokens) {
    if (grepl("^[0-9]+$", token)) output <- paste0(output, "[]")
    else output <- paste0(output, ".", decode_pointer(token))
  }
  output
}
attributes$schema_path <- vapply(attributes$json_pointer, pointer_to_schema,
                                 character(1))

# The original inventory's json_type reflects observed type(s), not a
# guaranteed type constraint. Multiple types per path must be retained.
previous <- read.csv(schema_file, stringsAsFactors = FALSE,
                     check.names = FALSE, na.strings = character())
if (!all(c("json_path", "json_type") %in% names(previous))) {
  stop("Schema CSV must contain json_path and json_type columns.")
}
current <- unique(attributes[c("schema_path", "node_type")])
names(current) <- c("json_path", "json_type")
current <- rbind(data.frame(json_path = "root", json_type = "array"), current)
current <- unique(current)
previous <- unique(previous[c("json_path", "json_type")])
previous$key <- paste(previous$json_path, previous$json_type, sep = "\t")
current$key <- paste(current$json_path, current$json_type, sep = "\t")

# When a field is consistently explicit null, the old schema's type may
# classify it differently. Flag discrepancies rather than discarding values.
# ------------------------------------------------------------
# Compare previous and current JSON schemas
# Safely handle zero-row differences
# ------------------------------------------------------------

only_previous <- previous[
  !(previous$key %in% current$key),
  c("json_path", "json_type"),
  drop = FALSE
]

only_current <- current[
  !(current$key %in% previous$key),
  c("json_path", "json_type"),
  drop = FALSE
]

schema_diff <- rbind(
  transform(
    only_previous,
    status = rep("in_previous_not_observed", nrow(only_previous))
  ),
  transform(
    only_current,
    status = rep("observed_not_in_previous", nrow(only_current))
  )
)

# Arrange columns consistently
schema_diff <- schema_diff[
  , c("status", "json_path", "json_type"),
  drop = FALSE
]

# Display validation results
cat("\nSchema comparison:\n")
cat("Fields in previous schema but not observed:",
    nrow(only_previous), "\n")
cat("Fields observed but absent from previous schema:",
    nrow(only_current), "\n")
if (!nrow(schema_diff)) schema_diff <- data.frame(
  status = character(), json_path = character(), json_type = character())
schema_diff <- schema_diff[order(schema_diff$status, schema_diff$json_path), ]

# Aggregate only: reports/validation should never contain patient-level values.
pattern_counts <- aggregate(attributes$patient_id,
                            by = list(json_path = attributes$schema_path,
                                      json_type = attributes$node_type),
                            FUN = length)
names(pattern_counts)[3] <- "total_nodes"
distinct_cases <- aggregate(attributes$patient_id,
                            by = list(json_path = attributes$schema_path,
                                      json_type = attributes$node_type),
                            FUN = function(v) length(unique(v)))
names(distinct_cases)[3] <- "cases_with_path"
field_catalogue <- merge(pattern_counts, distinct_cases,
                         by = c("json_path", "json_type"), sort = TRUE)

validation <- data.frame(
  metric = c("patient_rows", "unique_patient_ids", "coad", "read",
             "source_attribute_nodes", "observed_schema_path_type_pairs",
             "previous_schema_path_type_pairs", "schema_discrepancies",
             "source_copy_checksum_match", "observed_birth_year_count"),
  value = c(nrow(patient), length(unique(patient$patient_id)),
            unname(project_counts["TCGA-COAD"]), unname(project_counts["TCGA-READ"]),
            nrow(attributes), nrow(current), nrow(previous), nrow(schema_diff),
            1L, if ("demographic__year_of_birth" %in% names(patient))
              sum(!is.na(patient$demographic__year_of_birth)) else 0L),
  stringsAsFactors = FALSE)

# Write outputs only after validation. The existing .gitignore protects data/;
# confirm reports contain aggregate values only before committing to Git.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
copy_file <- file.path(out_dir, "clinical_source_exact_copy.json")
if (!file.copy(source_file, copy_file, overwrite = TRUE)) stop("Could not copy source JSON")
if (!identical(unname(tools::md5sum(source_file)),
               unname(tools::md5sum(copy_file)))) stop("JSON checksum mismatch")
write.csv(patient, file.path(out_dir, "patient.csv"), row.names = FALSE, na = "")
write.csv(attributes, file.path(out_dir, "clinical_source_attributes.csv"),
          row.names = FALSE, na = "")
write.csv(field_catalogue, file.path(report_dir, "01_clinical_field_catalogue.csv"),
          row.names = FALSE)
write.csv(schema_diff, file.path(report_dir, "01_schema_coverage_differences.csv"),
          row.names = FALSE)
write.csv(validation, file.path(report_dir, "01_patient_validation.csv"),
          row.names = FALSE)
cat("\nPhase 5 patient/source extraction completed.\n")
print(validation, row.names = FALSE)
cat("\nSchema discrepancies: ", nrow(schema_diff),
    " (inspect results/validation/01_schema_coverage_differences.csv)\n", sep = "")
cat("All nested clinical records are retained, but separate relational tables",
    " for diagnosis, treatment, follow-up etc. are not yet built.\n")
cat("OMOP PERSON has not yet been created.\n")
