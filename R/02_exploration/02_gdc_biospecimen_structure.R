# =============================================================================
# TCGA COAD/READ -> OMOP | GDC Release 46.0
# 02_gdc_biospecimen_structure.R
# Explore biospecimen TSV structure, identifiers and lineage.
# Inputs are local, unmodified GDC TSVs. Outputs contain aggregates only.
# Run from the tcga-crc-omop project root in RStudio.
# =============================================================================

input_dir <- "data/raw/gdc_r46/biospecimen"
output_dir <- "results/exploration"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

entities <- c("sample", "portion", "analyte", "aliquot", "slide")
id_cols <- c(sample = "samples.sample_id",
             portion = "portions.portion_id",
             analyte = "analytes.analyte_id",
             aliquot = "aliquots.aliquot_id",
             slide = "slides.slide_id")

# Read the source TSVs (read.delim preserves exact column names).
tables <- setNames(lapply(entities, function(entity) {
  path <- file.path(input_dir, paste0(entity, ".tsv"))
  if (!file.exists(path)) stop("Missing input file: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}), entities)

# Identify missing IDs (including blank values) without exporting identifiers.
missing_id <- function(x) is.na(x) | trimws(as.character(x)) == ""

# One row per entity? Slide export is expected to contain repeated slide IDs.
entity_summary <- do.call(rbind, lapply(entities, function(entity) {
  dat <- tables[[entity]]
  ids <- dat[[id_cols[[entity]]]]
  valid <- !missing_id(ids)
  data.frame(entity = entity,
             rows = nrow(dat),
             distinct_ids = length(unique(ids[valid])),
             missing_ids = sum(!valid),
             exact_duplicate_rows = sum(duplicated(dat)),
             stringsAsFactors = FALSE)
}))
rownames(entity_summary) <- NULL

cat("\n=== BIOSPECIMEN ENTITY SUMMARY ===\n")
print(entity_summary, row.names = FALSE)

# Check that child IDs are linked to exactly one distinct parent ID.
# This tests observed associations, not whether referenced parents exist.
parent_checks <- list(
  list(child = "portion", parent = "sample", parent_col = "samples.sample_id"),
  list(child = "analyte", parent = "portion", parent_col = "portions.portion_id"),
  list(child = "aliquot", parent = "analyte", parent_col = "analytes.analyte_id"),
  list(child = "slide", parent = "portion", parent_col = "portions.portion_id"),
  list(child = "slide", parent = "sample", parent_col = "samples.sample_id"),
  list(child = "sample", parent = "case", parent_col = "cases.case_id")
)

relationship_summary <- do.call(rbind, lapply(parent_checks, function(check) {
  dat <- tables[[check$child]]
  child_ids <- dat[[id_cols[[check$child]]]]
  parent_ids <- dat[[check$parent_col]]
  valid_child <- !missing_id(child_ids)
  valid_parent <- !missing_id(parent_ids)
  valid <- valid_child & valid_parent
  pairs <- unique(data.frame(child_id = child_ids[valid],
                             parent_id = parent_ids[valid],
                             stringsAsFactors = FALSE))
  n_parents <- if (nrow(pairs)) {
    as.integer(table(pairs$child_id))
  } else integer(0)
  data.frame(child = check$child, parent = check$parent,
             rows_missing_parent_id = sum(valid_child & !valid_parent),
             distinct_child_ids_with_parent = length(n_parents),
             child_ids_with_multiple_parents = sum(n_parents > 1L),
             max_parents_per_child = if (length(n_parents)) max(n_parents) else NA_integer_,
             duplicate_child_parent_pairs = sum(duplicated(data.frame(
               child_id = child_ids[valid], parent_id = parent_ids[valid]))),
             stringsAsFactors = FALSE)
}))
rownames(relationship_summary) <- NULL
cat("\n=== PARENT ASSOCIATION SUMMARY ===\n")
print(relationship_summary, row.names = FALSE)

# Referential integrity: every nonmissing foreign key should exist in its
# parent table. The clinical case IDs come from the clinical JSON, not here.
links <- list(
  list(child = "portion", parent = "sample", fk = "samples.sample_id"),
  list(child = "analyte", parent = "portion", fk = "portions.portion_id"),
  list(child = "aliquot", parent = "analyte", fk = "analytes.analyte_id"),
  list(child = "slide", parent = "portion", fk = "portions.portion_id"),
  list(child = "slide", parent = "sample", fk = "samples.sample_id")
)
integrity_summary <- do.call(rbind, lapply(links, function(link) {
  ids <- tables[[link$child]][[link$fk]]
  valid_ids <- ids[!missing_id(ids)]
  parent_ids <- tables[[link$parent]][[id_cols[[link$parent]]]]
  parent_ids <- parent_ids[!missing_id(parent_ids)]
  orphan <- unique(valid_ids[!valid_ids %in% parent_ids])
  data.frame(child = link$child, parent = link$parent,
             missing_fk_rows = sum(missing_id(ids)),
             unmatched_distinct_parent_ids = length(orphan),
             stringsAsFactors = FALSE)
}))
rownames(integrity_summary) <- NULL
cat("\n=== CROSS-TABLE REFERENTIAL INTEGRITY ===\n")
print(integrity_summary, row.names = FALSE)

# Slide-specific checks: repeated slide IDs, unique slide-portion pairs,
# distinct sample associations, and internal attribute consistency.
slide <- tables$slide
slide_ids <- slide[[id_cols[["slide"]]]]
slide_counts <- table(slide_ids[!missing_id(slide_ids)])
slide_frequency <- data.frame(
  occurrences_per_slide_id = as.integer(names(table(slide_counts))),
  number_of_slide_ids = as.integer(table(slide_counts)),
  stringsAsFactors = FALSE)
slide_cols <- grep("^slides\\.", names(slide), value = TRUE)
slide_attrs <- setdiff(slide_cols, "slides.slide_id")
repeated <- names(slide_counts)[slide_counts > 1]
repeated_rows <- slide[slide_ids %in% repeated, , drop = FALSE]
slide_attribute_conflicts <- if (length(repeated) == 0) 0L else {
  sum(vapply(split(repeated_rows, repeated_rows$`slides.slide_id`),
             function(group) any(vapply(slide_attrs, function(column) {
               length(unique(group[[column]])) > 1L
             }, logical(1))), logical(1)))
}
slide_pair <- slide[, c("slides.slide_id", "portions.portion_id"), drop = FALSE]
slide_summary <- data.frame(
  metric = c("slide_rows", "distinct_slide_ids", "repeated_slide_ids",
             "exact_duplicate_slide_rows", "duplicate_slide_portion_pairs",
             "repeated_slide_ids_with_attribute_conflicts"),
  value = c(nrow(slide), length(slide_counts), length(repeated),
            sum(duplicated(slide)), sum(duplicated(slide_pair)),
            slide_attribute_conflicts),
  stringsAsFactors = FALSE)
cat("\n=== SLIDE SUMMARY ===\n")
print(slide_summary, row.names = FALSE)
print(slide_frequency, row.names = FALSE)

# Save only aggregated results. No patient, sample or slide IDs are exported.
write.csv(entity_summary, file.path(output_dir, "02_biospecimen_entity_summary.csv"), row.names = FALSE)
write.csv(relationship_summary, file.path(output_dir, "02_biospecimen_relationship_summary.csv"), row.names = FALSE)
write.csv(integrity_summary, file.path(output_dir, "02_biospecimen_integrity_summary.csv"), row.names = FALSE)
write.csv(slide_summary, file.path(output_dir, "02_biospecimen_slide_summary.csv"), row.names = FALSE)
write.csv(slide_frequency, file.path(output_dir, "02_biospecimen_slide_frequency.csv"), row.names = FALSE)
cat("\nCompleted. Aggregate outputs saved in: ", output_dir, "\n", sep = "")



#  Clinical–biospecimen case linkage validation -------------------------

# Read GDC clinical JSON
clinical <- jsonlite::fromJSON(
  "data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json"
)

# Extract unique case IDs from clinical and biospecimen sources
clinical_ids <- unique(clinical$case_id)
biospecimen_ids <- unique(tables$sample$`cases.case_id`)

# Compare case identifiers between the two datasets
linkage_summary <- data.frame(
  metric = c(
    "Clinical cases",
    "Biospecimen cases",
    "Biospecimen cases absent from clinical",
    "Clinical cases without samples"
  ),
  value = c(
    length(clinical_ids),
    length(biospecimen_ids),
    length(setdiff(biospecimen_ids, clinical_ids)),
    length(setdiff(clinical_ids, biospecimen_ids))
  )
)

# Display results
cat("\n=== CLINICAL–BIOSPECIMEN CASE LINKAGE ===\n")
print(linkage_summary, row.names = FALSE)

# Save aggregate results
write.csv(
  linkage_summary,
  file.path(output_dir, "02_clinical_biospecimen_linkage.csv"),
  row.names = FALSE
)
