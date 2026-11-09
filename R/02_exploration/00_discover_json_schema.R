
# ==============================================================
# TCGA GDC - JSON Hierarchical Schema Discovery
#
# Objective:
# Discover all JSON field paths and their structural types.
#
# Does NOT:
# - Count patients or clinical events
# - Print individual patient values
# - Modify the original source file
# - Perform OMOP mapping
# ==============================================================

# STEP 1 - Load package ----------------------------------------

library(jsonlite)


# STEP 2 - Read the original JSON ------------------------------

json_path <- paste0(
  "data/raw/gdc_r46/clinical/",
  "clinical.cohort.2026-10-29.json"
)

stopifnot(file.exists(json_path))

# Preserve JSON objects and arrays as lists.
# Do not automatically simplify them into data frames.

raw_json <- fromJSON(
  json_path,
  simplifyVector = FALSE
)


# STEP 3 - Define a function to identify JSON types ------------

get_json_type <- function(x) {
  
  if (is.null(x)) {
    return("null")
  }
  
  if (is.list(x)) {
    
    # Unnamed lists represent JSON arrays.
    # Empty objects/arrays may require separate verification.
    
    if (is.null(names(x))) {
      return("array")
    }
    
    return("object")
  }
  
  if (is.character(x)) return("string")
  if (is.logical(x))   return("boolean")
  if (is.numeric(x))   return("number")
  
  return("other")
}


# STEP 4 - Discover every unique field path --------------------

# Store unique paths in an environment.
# Each path can have one or more observed JSON types.

schema <- new.env(
  hash = TRUE,
  parent = emptyenv()
)

add_schema_entry <- function(path, type) {
  
  if (exists(path, envir = schema, inherits = FALSE)) {
    
    old_types <- get(path, envir = schema)
    
    assign(
      path,
      unique(c(old_types, type)),
      envir = schema
    )
    
  } else {
    
    assign(path, type, envir = schema)
  }
}


walk_json <- function(x, path = "root") {
  
  type <- get_json_type(x)
  
  add_schema_entry(path, type)
  
  if (type == "object") {
    
    for (field in names(x)) {
      
      child_path <- paste0(path, ".", field)
      
      walk_json(x[[field]], child_path)
    }
    
  } else if (type == "array") {
    
    # [] means this level contains an array.
    # Inspect every element, not just the first.
    
    element_path <- paste0(path, "[]")
    
    for (element in x) {
      walk_json(element, element_path)
    }
  }
  
  invisible(NULL)
}


# STEP 5 - Run schema discovery --------------------------------

walk_json(raw_json)


# STEP 6 - Create the schema inventory -------------------------

paths <- sort(ls(schema))

schema_inventory <- data.frame(
  json_path = paths,
  
  json_type = vapply(
    paths,
    function(p) {
      paste(
        sort(get(p, envir = schema)),
        collapse = " | "
      )
    },
    character(1)
  ),
  
  stringsAsFactors = FALSE
)


# STEP 7 - Display the hierarchical structure -----------------

# Indentation reflects the nesting depth.
# [] remains visible to identify arrays.

cat("\n=== GDC JSON HIERARCHY ===\n\n")

tree_lines <- vapply(
  seq_len(nrow(schema_inventory)),
  function(i) {
    
    path <- schema_inventory$json_path[i]
    type <- schema_inventory$json_type[i]
    
    parts <- strsplit(path, ".", fixed = TRUE)[[1]]
    
    depth <- length(parts) - 1
    
    label <- tail(parts, 1)
    
    paste0(
      strrep("  ", depth),
      label,
      "  [", type, "]"
    )
  },
  character(1)
)

cat(paste(tree_lines, collapse = "\n"), "\n")


# STEP 8 - Save reports ----------------------------------------

output_dir <- "results/exploration"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  schema_inventory,
  file.path(output_dir, "01_gdc_json_schema.csv"),
  row.names = FALSE,
  na = ""
)

writeLines(
  c("GDC JSON HIERARCHY", "", tree_lines),
  file.path(output_dir, "01_gdc_json_schema_tree.txt")
)


# STEP 9 - Completion message ----------------------------------

cat("\nSchema discovery complete.\n")

cat(
  "Unique structural paths:",
  nrow(schema_inventory),
  "\n"
)

cat("\nReports saved to:\n")
cat("  results/exploration/01_gdc_json_schema.csv\n")
cat("  results/exploration/01_gdc_json_schema_tree.txt\n")
