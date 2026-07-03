# Restore the uploaded ASV abundance matrix from the compact archive in data/.
#
# Archive format (ASVC1):
#   5-byte magic string "ASVC1"
#   unsigned LEB128 varint: number of ASVs
#   unsigned LEB128 varint: number of samples
#   ASV numeric identifiers stored as successive varint differences
#   for each sample: number of non-zero entries, then row-index differences
#   and integer counts, all encoded as unsigned varints
#
# The binary stream is compressed with XZ and split into four Base64 text files
# so it can be stored reliably through GitHub's text-file interface.

restore_abundance_from_archive <- function(output_path, data_dir) {
  parts <- sort(Sys.glob(file.path(data_dir, "asv_counts.compact.b64.part*")))
  if (length(parts) != 4L) {
    stop(
      "Expected four compact abundance archive parts in data/, but found ",
      length(parts), ".",
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
  archive_raw <- tryCatch(
    memDecompress(compressed_raw, type = "xz"),
    error = function(e) {
      stop("Unable to decompress the ASV abundance archive: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (length(archive_raw) < 8L || rawToChar(archive_raw[1:5]) != "ASVC1") {
    stop("The abundance archive has an invalid ASVC1 header.", call. = FALSE)
  }

  position <- 6L
  read_varint <- function() {
    value <- 0
    shift <- 0
    repeat {
      if (position > length(archive_raw)) {
        stop("Unexpected end of the abundance archive.", call. = FALSE)
      }
      byte <- as.integer(archive_raw[position])
      position <<- position + 1L
      value <- value + bitwAnd(byte, 127L) * (2 ^ shift)
      if (bitwAnd(byte, 128L) == 0L) break
      shift <- shift + 7L
      if (shift > 49L) stop("Invalid varint in the abundance archive.", call. = FALSE)
    }
    as.integer(value)
  }

  n_asv <- read_varint()
  n_samples <- read_varint()
  expected_samples <- c(
    "RSV1", "RSV2", "RSV3", "RRV1", "RRV2", "RRV3",
    "SSV1", "SSV2", "SSV3", "SRV1", "SRV2", "SRV3",
    "ESV1", "ESV2", "ESV3", "ERV1", "ERV2", "ERV3"
  )

  if (n_asv != 24122L || n_samples != length(expected_samples)) {
    stop(
      "Unexpected archive dimensions: ", n_asv, " ASVs x ", n_samples,
      " samples.",
      call. = FALSE
    )
  }

  asv_numbers <- integer(n_asv)
  previous_id <- 0L
  for (i in seq_len(n_asv)) {
    previous_id <- previous_id + read_varint()
    asv_numbers[i] <- previous_id
  }

  counts <- matrix(0L, nrow = n_asv, ncol = n_samples)
  colnames(counts) <- expected_samples

  for (j in seq_len(n_samples)) {
    nonzero_count <- read_varint()
    previous_row <- 0L
    if (nonzero_count > 0L) {
      for (k in seq_len(nonzero_count)) {
        row_zero_based <- previous_row + read_varint()
        count_value <- read_varint()
        row_one_based <- row_zero_based + 1L
        if (row_one_based < 1L || row_one_based > n_asv || count_value < 0L) {
          stop("Invalid row index or count in the abundance archive.", call. = FALSE)
        }
        counts[row_one_based, j] <- count_value
        previous_row <- row_zero_based
      }
    }
  }

  if (position != length(archive_raw) + 1L) {
    stop("Unexpected trailing bytes in the abundance archive.", call. = FALSE)
  }

  expected_depth <- 31556L
  observed_depth <- colSums(counts)
  if (any(observed_depth != expected_depth)) {
    bad <- expected_samples[observed_depth != expected_depth]
    stop(
      "Restored sequencing depth is incorrect for: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  output <- data.frame(
    `#ID` = paste("ASV", asv_numbers),
    as.data.frame(counts, check.names = FALSE),
    check.names = FALSE
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    output,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    fileEncoding = "UTF-8"
  )

  message(
    "Restored ASV abundance matrix: ", n_asv, " ASVs x ", n_samples,
    " samples -> ", output_path
  )
  invisible(output_path)
}

if (exists("PATHS") && !file.exists(PATHS$abundance)) {
  restore_abundance_from_archive(PATHS$abundance, DATA_DIR)
}
