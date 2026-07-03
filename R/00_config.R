# ---------------------------
# Global configuration
# ---------------------------

PROJECT_DIR <- normalizePath(".", winslash = "/", mustWork = FALSE)
DATA_DIR    <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR  <- file.path(PROJECT_DIR, "output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PATHS <- list(
  metadata  = file.path(DATA_DIR, "sample_metadata.tsv"),
  taxonomy  = file.path(DATA_DIR, "taxonomy.tsv"),
  abundance = file.path(DATA_DIR, "asv_abundance.tsv")
)

# The uploaded sample codes occur as paired groups: RSV/RRV, SSV/SRV, ESV/ERV.
# This project interprets the second character as compartment:
# S = Rhizosphere, R = Root. Change sample_metadata.tsv if this is not correct.
TREATMENT_LEVELS <- c("R", "S", "E")
TREATMENT_LABELS <- c(R = "R", S = "S", E = "E")

# Palette chosen to match the reference figure: ochre, blue, purple.
TREATMENT_COLORS <- c(
  R = "#D6A23B",
  S = "#56A6C9",
  E = "#6E5195"
)

# Differential-abundance comparisons used in panels d/e.
# Each row is group1 vs group2; the displayed difference is group1 - group2.
COMPARISONS <- data.frame(
  comparison_id = c("R_vs_S", "R_vs_E"),
  group1 = c("R", "R"),
  group2 = c("S", "E"),
  stringsAsFactors = FALSE
)

TOP_N_FAMILIES <- 20
P_ADJUST_METHOD <- "BH"
PERMUTATIONS <- 999
RANDOM_SEED <- 20260703

FIG_WIDTHS <- list(
  panel_a = 8.4,
  panel_b = 9.2,
  panel_c = 8.6,
  panel_de = 10.5,
  combined = 13.5
)
