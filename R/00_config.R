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

# Reconstruct the uploaded abundance and family-level taxonomy tables
# automatically when only the compact archive parts are present.
source(file.path(PROJECT_DIR, "R", "00_restore_abundance.R"))
source(file.path(PROJECT_DIR, "R", "00_restore_taxonomy.R"))

# The uploaded sample codes occur as paired groups: RSV/RRV, SSV/SRV, ESV/ERV.
# This project interprets the second character as compartment:
# S = Rhizosphere, R = Root. Change sample_metadata.tsv if this is not correct.
TREATMENT_LEVELS <- c("R", "S", "E")
TREATMENT_LABELS <- c(R = "R", S = "S", E = "E")

# Palette chosen to match the supplied reference figures.
TREATMENT_COLORS <- c(
  R = "#D6A23B",
  S = "#56A6C9",
  E = "#6E5195"
)

# Differential-abundance comparisons used in Figure 2 panels d/e.
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

# ---------------------------
# Figure 3: network robustness and NetShift-style analysis
# ---------------------------

# R is treated as the focal/case group and S or E as the reference/control.
# Reverse case/control here if the biological treatment meanings are different.
NETSHIFT_COMPARISONS <- data.frame(
  comparison_id = c("R_vs_S", "R_vs_E"),
  case = c("R", "R"),
  control = c("S", "E"),
  case_label = c("R", "R"),
  control_label = c("S", "E"),
  stringsAsFactors = FALSE
)

NETWORK_TAXON_LEVEL <- "family"
NETWORK_TOP_N <- 50
NETWORK_MIN_PREVALENCE <- 0.30
NETWORK_PSEUDOCOUNT <- 0.5

# With the current n = 3 samples per treatment-compartment group, ordinary
# correlations are singular. The workflow therefore uses CLR transformation
# followed by shrinkage correlation. Increase replication before publication.
NETWORK_EDGE_QUANTILE <- 0.82
NETWORK_MIN_ABS_COR <- 0.20
NETWORK_MIN_EDGES <- 25
NETWORK_MAX_EDGES <- 140
NETWORK_POSITIVE_ONLY <- TRUE

ROBUSTNESS_REMOVAL_PROPORTIONS <- seq(0, 0.9, by = 0.1)
ROBUSTNESS_ITERATIONS <- 100

# NetShift-style driver rule: positive change in scaled betweenness plus a
# directional Neighbor Shift score above the stated threshold.
NESH_DRIVER_MIN <- 1.0
NESH_FALLBACK_TOP_N <- 3

NETSHIFT_EDGE_COLORS <- c(
  "Control only" = "#D73027",
  "Case only" = "#1A9850",
  "Both" = "#4575B4"
)

UPSET_MAX_COMBINATIONS <- 15

FIG3_WIDTHS <- list(
  panel_a = 5.2,
  panel_b = 11.2,
  panel_c = 7.0,
  combined = 15.5
)
