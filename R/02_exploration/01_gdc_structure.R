
# ==============================================================================
# 01_gdc_structure.R
#
# Project : TCGA COAD/READ to OMOP
# Purpose : Explore the structure of GDC clinical data (Release 46.0)
# Inputs  : GDC clinical JSON and clinical.tsv
# Outputs : Aggregate CSV summaries in results/exploration/
#
# No patient-level clinical data are exported.
# ==============================================================================

library(jsonlite)
library(dplyr)
library(readr)

# 1. File paths ---------------------------------------------------------------

json_path <- "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
tsv_path  <- "data/raw/gdc_r46/clinical/clinical.tsv"
out_dir   <- "results/exploration"

stopifnot(file.exists(json_path), file.exists(tsv_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Read clinical data -------------------------------------------------------

gdc <- fromJSON(json_path)
clinical <- read.delim(tsv_path, check.names = FALSE)

cat("\n=== CLINICAL DATA STRUCTURE ===\n")

n_cases <- length(gdc$submitter_id)
n_tsv_rows <- nrow(clinical)
n_tsv_columns <- ncol(clinical)
n_unique_tsv_cases <- length(unique(clinical$cases.submitter_id))

cat("JSON cases:", n_cases, "\n")
cat("TSV rows:", n_tsv_rows, "\n")
cat("TSV columns:", n_tsv_columns, "\n")
cat("Unique TSV cases:", n_unique_tsv_cases, "\n")

stopifnot(n_cases == n_unique_tsv_cases)

# 3. Repeated clinical TSV rows ----------------------------------------------

rows_per_patient <- table(clinical$cases.submitter_id)

tsv_distribution <- as.data.frame(
  table(as.integer(rows_per_patient))
)

names(tsv_distribution) <- c("rows_per_patient", "n_patients")
tsv_distribution$rows_per_patient <- as.integer(
  as.character(tsv_distribution$rows_per_patient)
)

cat("\n=== TSV ROWS PER PATIENT ===\n")
print(tsv_distribution)

# 4. Diagnosis counts ---------------------------------------------------------

diagnosis_counts <- vapply(
  gdc$diagnoses,
  function(x) if (is.null(x)) 0L else as.integer(NROW(x)),
  integer(1)
)

diagnosis_distribution <- as.data.frame(table(diagnosis_counts))
names(diagnosis_distribution) <- c("n_diagnoses", "n_patients")

cat("\n=== DIAGNOSES PER PATIENT ===\n")
print(diagnosis_distribution)

# 5. Diagnosis classifications ------------------------------------------------

diagnosis_types <- unlist(
  lapply(gdc$diagnoses, function(x) x$classification_of_tumor),
  use.names = FALSE
)

diagnosis_classification <- as.data.frame(
  table(diagnosis_types, useNA = "ifany")
)

names(diagnosis_classification) <- c("classification", "n_records")

cat("\n=== DIAGNOSIS CLASSIFICATIONS ===\n")
print(diagnosis_classification)

# 6. Primary diagnosis counts -------------------------------------------------

primary_counts <- vapply(
  gdc$diagnoses,
  function(x) {
    if (is.null(x)) return(0L)
    as.integer(sum(x$classification_of_tumor == "primary",
                   na.rm = TRUE))
  },
  integer(1)
)

primary_distribution <- as.data.frame(table(primary_counts))
names(primary_distribution) <- c("n_primary_diagnoses", "n_patients")

cat("\n=== PRIMARY DIAGNOSES PER PATIENT ===\n")
print(primary_distribution)

cat("Cases without a primary label:", sum(primary_counts == 0), "\n")

# 7. Treatment records --------------------------------------------------------

# Treatment tables have different columns.
# bind_rows() accommodates these differences by filling missing fields with NA.

treatment_tables <- unlist(
  lapply(gdc$diagnoses, function(d) d$treatments),
  recursive = FALSE
)

all_treatments <- bind_rows(treatment_tables)

treatment_counts <- vapply(
  gdc$diagnoses,
  function(diagnoses) {
    if (is.null(diagnoses) ||
        !"treatments" %in% names(diagnoses)) return(0L)
    
    as.integer(sum(vapply(
      diagnoses$treatments,
      function(t) if (is.null(t)) 0L else as.integer(NROW(t)),
      integer(1)
    )))
  },
  integer(1)
)

cat("\n=== TREATMENT RECORDS PER PATIENT ===\n")
print(summary(treatment_counts))

treatment_indicators <- as.data.frame(
  table(all_treatments$treatment_or_therapy, useNA = "ifany")
)

names(treatment_indicators) <- c("treatment_indicator", "n_records")

cat("\n=== TREATMENT INDICATORS ===\n")
print(treatment_indicators)

# Patients with at least one treatment record marked yes

patients_with_treatment_yes <- vapply(
  gdc$diagnoses,
  function(diagnoses) {
    if (is.null(diagnoses)) return(FALSE)
    
    any(vapply(
      diagnoses$treatments,
      function(t) {
        !is.null(t) &&
          any(t$treatment_or_therapy == "yes", na.rm = TRUE)
      },
      logical(1)
    ))
  },
  logical(1)
)

cat("Patients with treatment marked yes:",
    sum(patients_with_treatment_yes), "\n")

cat("Patients with no treatment marked yes:",
    sum(!patients_with_treatment_yes), "\n")

# 8. Treatment types marked yes -----------------------------------------------

positive_treatments <- all_treatments |>
  filter(treatment_or_therapy == "yes")

treatment_types <- positive_treatments |>
  count(treatment_type, name = "n_records", sort = TRUE)

cat("\n=== POSITIVE TREATMENT TYPES ===\n")
print(treatment_types, n = Inf)

# 9. Chemotherapy agent availability ------------------------------------------

chemo <- positive_treatments |>
  filter(treatment_type == "Chemotherapy")

agent_text <- trimws(as.character(chemo$therapeutic_agents))

agent_missing <- is.na(agent_text) |
  tolower(agent_text) %in%
  c("", "unknown", "not reported", "not available", "--")

chemo_agent_summary <- data.frame(
  category = c("Agent recorded", "Agent missing or unreported"),
  n_records = c(sum(!agent_missing), sum(agent_missing))
)

cat("\n=== CHEMOTHERAPY AGENT AVAILABILITY ===\n")
print(chemo_agent_summary)

# This is a preliminary completeness check.
# It does not validate drug names or OMOP vocabulary mappings.

# 10. Follow-up records -------------------------------------------------------

followup_counts <- vapply(
  gdc$follow_ups,
  function(x) if (is.null(x)) 0L else as.integer(NROW(x)),
  integer(1)
)

cat("\n=== FOLLOW-UP RECORDS PER PATIENT ===\n")
print(summary(followup_counts))

cat("Patients without follow-up records:",
    sum(followup_counts == 0), "\n")

# 11. Summary metrics ---------------------------------------------------------

summary_metrics <- data.frame(
  metric = c(
    "JSON cases",
    "TSV rows",
    "TSV columns",
    "Unique TSV cases",
    "Total diagnoses",
    "Patients with multiple diagnoses",
    "Patients without primary diagnosis label",
    "Total treatment records",
    "Treatment records marked yes",
    "Patients with treatment marked yes",
    "Patients without treatment records",
    "Patients without follow-up records",
    "Maximum diagnoses per patient",
    "Maximum treatment records per patient",
    "Maximum follow-up records per patient"
  ),
  value = c(
    n_cases,
    n_tsv_rows,
    n_tsv_columns,
    n_unique_tsv_cases,
    sum(diagnosis_counts),
    sum(diagnosis_counts > 1),
    sum(primary_counts == 0),
    nrow(all_treatments),
    nrow(positive_treatments),
    sum(patients_with_treatment_yes),
    sum(treatment_counts == 0),
    sum(followup_counts == 0),
    max(diagnosis_counts),
    max(treatment_counts),
    max(followup_counts)
  )
)

cat("\n=== SUMMARY ===\n")
print(summary_metrics)

# 12. Save aggregate results --------------------------------------------------

write_csv(summary_metrics,
          file.path(out_dir, "01_gdc_structure_summary.csv"))

write_csv(tsv_distribution,
          file.path(out_dir, "01_tsv_row_distribution.csv"))

write_csv(diagnosis_distribution,
          file.path(out_dir, "01_diagnosis_distribution.csv"))

write_csv(diagnosis_classification,
          file.path(out_dir, "01_diagnosis_classifications.csv"))

write_csv(primary_distribution,
          file.path(out_dir, "01_primary_diagnosis_distribution.csv"))

write_csv(treatment_indicators,
          file.path(out_dir, "01_treatment_indicators.csv"))

write_csv(treatment_types,
          file.path(out_dir, "01_positive_treatment_types.csv"))

write_csv(chemo_agent_summary,
          file.path(out_dir, "01_chemo_agent_availability.csv"))

cat("\nGDC clinical structure exploration completed.\n")
