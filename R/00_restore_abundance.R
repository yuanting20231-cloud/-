# Restore the uploaded abundance matrix from the compact archive stored in GitHub.
# The original matrix is 24,122 ASVs x 18 samples. It is stored as two
# Base64-encoded XZ parts to stay within GitHub connector text-file limits.

restore_abundance_from_archive <- function(output_path, data_dir) {
  parts <- sort(Sys.glob(file.path(data_dir, "asv_abundance.tsv.xz.b64.part*")))
  if (length(parts) == 0) {
    stop(
      "Neither data/asv_abundance.tsv nor archive parts were found.",
      call. = FALSE
    )
  }

  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop(
      "Package 'base64enc' is required to restore the abundance matrix. ",
      "Run source('install_packages.R') first.",
      call. = FALSE
    )
  }

  encoded <- paste0(
    vapply(
      parts,
      function(path) paste0(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = ""),
      character(1)
    ),
    collapse = ""
  )

  compressed_raw <- base64enc::base64decode(encoded)
  table_raw <- memDecompress(compressed_raw, type = "xz")

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  con <- file(output_path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(table_raw, con)

  if (file.info(output_path)$size != 1107200) {
    stop(
      "Restored abundance matrix has an unexpected file size; archive may be incomplete.",
      call. = FALSE
    )
  }

  invisible(output_path)
}

if (exists("PATHS") && !file.exists(PATHS$abundance)) {
  restore_abundance_from_archive(PATHS$abundance, DATA_DIR)
}
