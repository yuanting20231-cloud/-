packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "patchwork", "vegan", "FSA", "multcompView"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, dependencies = TRUE)
