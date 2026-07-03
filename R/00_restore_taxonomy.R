# Restore the uploaded family-level taxonomy table from the compact archive.
#
# Archive format (TAXF1):
#   5-byte magic string "TAXF1"
#   unsigned LEB128 varint: number of ASVs
#   unsigned LEB128 varint: number of unique family labels
#   family dictionary: UTF-8 byte length plus bytes for each label
#   for each ASV: numeric ASV-ID difference, then zero-based family index
#
# The binary stream is compressed with XZ and Base64 encoded across four files.

restore_taxonomy_from_archive <- function(output_path, data_dir) {
  parts <- sort(Sys.glob(file.path(data_dir, "taxonomy_family.compact.b64.part*")))
  if (length(parts) != 4L) {
    stop(
      "Expected four family-taxonomy archive parts in data/, but found ",
      length(parts), ".",
      call. = FALSE
    )
  }

  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop(
      "Package 'base64enc' is required to restore taxonomy.tsv. ",
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
      stop("Unable to decompress the family-taxonomy archive: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (length(archive_raw) < 8L || rawToChar(archive_raw[1:5]) != "TAXF1") {
    stop("The family-taxonomy archive has an invalid TAXF1 header.", call. = FALSE)
  }

  position <- 6L
  read_varint <- function() {
    value <- 0
    shift <- 0
    repeat {
      if (position > length(archive_raw)) {
        stop("Unexpected end of the family-taxonomy archive.", call. = FALSE)
      }
      byte <- as.integer(archive_raw[position])
      position <<- position + 1L
      value <- value + bitwAnd(byte, 127L) * (2 ^ shift)
      if (bitwAnd(byte, 128L) == 0L) break
      shift <- shift + 7L
      if (shift > 49L) stop("Invalid varint in the family-taxonomy archive.", call. = FALSE)
    }
    as.integer(value)
  }

  read_utf8_string <- function() {
    n_bytes <- read_varint()
    if (n_bytes == 0L) return("")
    end_position <- position + n_bytes - 1L
    if (end_position > length(archive_raw)) {
      stop("Unexpected end while reading a taxonomy label.", call. = FALSE)
    }
    value <- rawToChar(archive_raw[position:end_position])
    position <<- end_position + 1L
    enc2utf8(value)
  }

  n_asv <- read_varint()
  n_families <- read_varint()
  if (n_asv != 24122L || n_families < 1L) {
    stop(
      "Unexpected family-taxonomy dimensions: ", n_asv,
      " ASVs and ", n_families, " families.",
      call. = FALSE
    )
  }

  family_dictionary <- character(n_families)
  for (i in seq_len(n_families)) {
    family_dictionary[i] <- read_utf8_string()
  }

  asv_numbers <- integer(n_asv)
  family_index <- integer(n_asv)
  previous_id <- 0L
  for (i in seq_len(n_asv)) {
    previous_id <- previous_id + read_varint()
    asv_numbers[i] <- previous_id
    family_index[i] <- read_varint() + 1L
  }

  if (position != length(archive_raw) + 1L) {
    stop("Unexpected trailing bytes in the family-taxonomy archive.", call. = FALSE)
  }
  if (any(family_index < 1L | family_index > n_families)) {
    stop("Invalid family dictionary index in the taxonomy archive.", call. = FALSE)
  }

  output <- data.frame(
    `#ID` = paste("ASV", asv_numbers),
    family = family_dictionary[family_index],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (anyDuplicated(output$`#ID`) || anyNA(output$family)) {
    stop("Restored taxonomy contains duplicate ASV IDs or missing family labels.", call. = FALSE)
  }

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
    "Restored family-level taxonomy: ", n_asv, " ASVs, ", n_families,
    " family labels -> ", output_path
  )
  invisible(output_path)
}

if (exists("PATHS") && !file.exists(PATHS$taxonomy)) {
  restore_taxonomy_from_archive(PATHS$taxonomy, DATA_DIR)
}
