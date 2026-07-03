packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "patchwork", "vegan", "FSA", "multcompView", "base64enc",
  "igraph", "tidygraph", "ggraph", "corpcor", "scales"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, dependencies = TRUE)
