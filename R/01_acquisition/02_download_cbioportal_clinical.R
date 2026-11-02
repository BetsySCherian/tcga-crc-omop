# ==============================================================================
# 02_download_cbioportal_clinical.R
# Purpose : Download cBioPortal TCGA PanCanAtlas COAD/READ clinical tables,
#           pinned to a fixed cBioPortal datahub commit, and verify checksums.
# Role    : Source reconciliation only. GDC Release 46.0 remains the ETL source.
# Needs   : install.packages("openssl")
# ==============================================================================

study_id       <- "coadread_tcga_pan_can_atlas_2018"
datahub_commit <- "e1b75cd25bf7436d16719b7d76117ad8f3f66e54"   # cBioPortal/datahub, 2026-10-27
out_dir        <- file.path("data", "raw", "cbioportal", study_id)

# Expected file fingerprints, taken from the datahub at the pinned commit
expected <- data.frame(
  file   = c("data_clinical_patient.txt", "data_clinical_sample.txt"),
  bytes  = c(175194, 105858),
  sha256 = c("3c5738f41b9047313af8f809c244bc66419991379b0c1978816dfc56c184c313",
             "b9cbcc7cbc7471b44a3767aac48d4d0625001e3c07f8776c50f77364341f7817")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- sprintf(
  "https://media.githubusercontent.com/media/cBioPortal/datahub/%s/public/%s",
  datahub_commit, study_id
)

for (i in seq_len(nrow(expected))) {
  f    <- expected$file[i]
  dest <- file.path(out_dir, f)

  download.file(paste0(base_url, "/", f), dest, mode = "wb", quiet = TRUE)

  got_bytes  <- file.size(dest)
  got_sha256 <- as.character(openssl::sha256(file(dest)))

  if (got_bytes != expected$bytes[i])    stop("Size mismatch for ", f, ": ", got_bytes, " bytes")
  if (got_sha256 != expected$sha256[i])  stop("Checksum mismatch for ", f)

  message("OK  ", f, "  (", got_bytes, " bytes, SHA-256 verified)")
}
