`%||%` <- function(x, y) if (is.null(x)) y else x

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "patchwork", "vegan", "FSA", "multcompView"
)

check_packages <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nInstall them with: install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
}

check_packages()
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(FSA)
  library(multcompView)
})

read_utf_table <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  enc <- readr::guess_encoding(path, n_max = 1000)$encoding[1] %||% "UTF-8"
  suppressMessages(readr::read_tsv(path, locale = locale(encoding = enc), show_col_types = FALSE, progress = FALSE))
}

clean_taxon <- function(x, fallback = "Unclassified") {
  x <- as.character(x)
  bad <- is.na(x) | x == "" | grepl("^(unclassified|uncultured)( |$)", x, ignore.case = TRUE)
  x[bad] <- fallback
  x
}

load_microbiome_data <- function(paths = PATHS) {
  check_packages()

  metadata <- read_utf_table(paths$metadata) %>%
    mutate(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Treatment = factor(as.character(Treatment), levels = TREATMENT_LEVELS),
      Compartment = factor(as.character(Compartment), levels = c("Rhizosphere", "Root"))
    )

  taxonomy <- read_utf_table(paths$taxonomy)
  if (!"#ID" %in% names(taxonomy)) names(taxonomy)[1] <- "#ID"
  taxonomy <- taxonomy %>%
    mutate(`#ID` = str_squish(as.character(`#ID`))) %>%
    distinct(`#ID`, .keep_all = TRUE)

  abundance_raw <- read_utf_table(paths$abundance)
  if (!"#ID" %in% names(abundance_raw)) names(abundance_raw)[1] <- "#ID"
  abundance_raw <- abundance_raw %>%
    mutate(`#ID` = str_squish(as.character(`#ID`))) %>%
    distinct(`#ID`, .keep_all = TRUE)

  sample_cols <- intersect(metadata$Sample, names(abundance_raw))
  if (length(sample_cols) != nrow(metadata)) {
    missing_samples <- setdiff(metadata$Sample, names(abundance_raw))
    stop(
      "Abundance table is missing sample columns: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  abundance <- abundance_raw %>%
    select(`#ID`, all_of(metadata$Sample)) %>%
    mutate(across(-`#ID`, ~ suppressWarnings(as.numeric(.x))))

  if (anyNA(abundance[-1])) {
    stop("Non-numeric or missing values were found in asv_abundance.tsv.", call. = FALSE)
  }
  if (any(as.matrix(abundance[-1]) < 0)) {
    stop("Negative abundance values are not allowed.", call. = FALSE)
  }

  shared_ids <- intersect(abundance$`#ID`, taxonomy$`#ID`)
  if (length(shared_ids) == 0) {
    stop("No matching ASV IDs between abundance and taxonomy tables.", call. = FALSE)
  }

  abundance <- abundance %>% filter(`#ID` %in% shared_ids)
  taxonomy <- taxonomy %>% filter(`#ID` %in% shared_ids)

  sample_taxa <- abundance %>%
    column_to_rownames("#ID") %>%
    as.matrix() %>%
    t()
  storage.mode(sample_taxa) <- "numeric"
  sample_taxa <- sample_taxa[metadata$Sample, , drop = FALSE]

  if (any(rowSums(sample_taxa) == 0)) {
    stop("At least one sample has a total abundance of zero.", call. = FALSE)
  }

  list(
    metadata = metadata,
    taxonomy = taxonomy,
    abundance = abundance,
    sample_taxa = sample_taxa
  )
}

relative_abundance_matrix <- function(sample_taxa) {
  sweep(sample_taxa, 1, rowSums(sample_taxa), "/")
}

theme_reference <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "plain", hjust = 0.5),
      legend.title = element_blank(),
      legend.key.height = grid::unit(0.38, "cm"),
      panel.spacing = grid::unit(0.7, "lines")
    )
}

save_plot <- function(plot, filename, width, height, dpi = 600) {
  ggsave(
    filename = file.path(OUTPUT_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

dunn_letters <- function(df, response, group = "Treatment") {
  response_vec <- df[[response]]
  group_vec <- droplevels(factor(df[[group]]))
  out <- setNames(rep("a", nlevels(group_vec)), levels(group_vec))

  if (nlevels(group_vec) < 2 || length(unique(response_vec)) < 2) return(out)

  dt <- tryCatch(
    FSA::dunnTest(response_vec, group_vec, method = "bh")$res,
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(out)

  pvals <- dt$P.adj
  names(pvals) <- gsub(" - ", "-", dt$Comparison, fixed = TRUE)
  letters <- multcompView::multcompLetters(pvals, threshold = 0.05)$Letters
  out[names(letters)] <- letters
  out
}

family_relative_long <- function(dat) {
  taxonomy_family <- dat$taxonomy %>%
    transmute(
      `#ID`,
      Family = clean_taxon(family, fallback = "Unclassified")
    )

  dat$abundance %>%
    pivot_longer(-`#ID`, names_to = "Sample", values_to = "Count") %>%
    left_join(taxonomy_family, by = "#ID") %>%
    group_by(Sample, Family) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(RelAbund = 100 * Count / sum(Count)) %>%
    ungroup() %>%
    left_join(dat$metadata, by = "Sample")
}

top_families <- function(family_long, n = TOP_N_FAMILIES) {
  family_long %>%
    group_by(Family) %>%
    summarise(Mean = mean(RelAbund), .groups = "drop") %>%
    arrange(desc(Mean)) %>%
    slice_head(n = n) %>%
    pull(Family)
}

student_difference <- function(x, y, conf.level = 0.95) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2 || n2 < 2) {
    return(tibble(diff = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  }

  test <- t.test(x, y, var.equal = TRUE, conf.level = conf.level)
  tibble(
    diff = mean(x) - mean(y),
    lower = unname(test$conf.int[1]),
    upper = unname(test$conf.int[2]),
    p = test$p.value
  )
}

prepare_family_comparison <- function(family_long, compartment, group1, group2, top_n = TOP_N_FAMILIES) {
  sub <- family_long %>%
    filter(Compartment == compartment, Treatment %in% c(group1, group2))

  top <- top_families(sub, n = top_n)
  sub <- sub %>% filter(Family %in% top)

  stats <- lapply(top, function(fam) {
    d <- sub %>% filter(Family == fam)
    x <- d %>% filter(Treatment == group1) %>% pull(RelAbund)
    y <- d %>% filter(Treatment == group2) %>% pull(RelAbund)
    student_difference(x, y) %>% mutate(Family = fam)
  }) %>% bind_rows() %>%
    mutate(p_adj = p.adjust(p, method = P_ADJUST_METHOD))

  means <- sub %>%
    group_by(Family, Treatment) %>%
    summarise(Mean = mean(RelAbund), SD = sd(RelAbund), .groups = "drop")

  result <- stats %>%
    left_join(
      means %>% select(Family, Treatment, Mean) %>%
        pivot_wider(names_from = Treatment, values_from = Mean, names_prefix = "mean_"),
      by = "Family"
    ) %>%
    arrange(p, desc(abs(diff)))

  significant <- result %>% filter(!is.na(p), p < 0.05)
  if (nrow(significant) == 0) {
    warning(
      "No unadjusted P < 0.05 families for ", compartment, " ", group1, " vs ", group2,
      "; displaying the five smallest P values."
    )
    significant <- result %>% slice_head(n = min(5, nrow(result)))
  }

  order_levels <- rev(significant$Family)
  list(
    stats = significant %>% mutate(Family = factor(Family, levels = order_levels)),
    means = means %>%
      filter(Family %in% significant$Family) %>%
      mutate(Family = factor(Family, levels = order_levels)),
    group1 = group1,
    group2 = group2,
    compartment = compartment
  )
}

plot_family_comparison <- function(prepared) {
  g1 <- prepared$group1
  g2 <- prepared$group2
  stats <- prepared$stats
  means <- prepared$means

  p_mean <- ggplot(means, aes(x = Mean, y = Family, fill = Treatment)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.68) +
    scale_fill_manual(
      values = TREATMENT_COLORS[c(g1, g2)],
      breaks = c(g1, g2),
      labels = TREATMENT_LABELS[c(g1, g2)]
    ) +
    labs(x = "Mean proportion (%)", y = NULL) +
    theme_reference(9) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(size = 9)
    )

  x_range <- range(c(stats$lower, stats$upper, 0), na.rm = TRUE)
  if (!all(is.finite(x_range)) || diff(x_range) == 0) x_range <- c(-1, 1)
  p_diff <- ggplot(stats, aes(x = diff, y = Family)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey40") +
    geom_segment(aes(x = lower, xend = upper, y = Family, yend = Family),
                 linewidth = 0.45, colour = "grey45") +
    geom_point(aes(colour = diff > 0), size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = TREATMENT_COLORS[[g1]], `FALSE` = TREATMENT_COLORS[[g2]]), guide = "none") +
    scale_x_continuous(limits = x_range, expand = expansion(mult = c(0.08, 0.12))) +
    labs(
      x = paste0("Difference in mean proportions (%)\n", g1, " - ", g2),
      y = NULL,
      title = "95% confidence intervals"
    ) +
    theme_reference(9) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 9)
    )

  p_p <- ggplot(stats, aes(x = 1, y = Family, label = sprintf("%.4f", p))) +
    geom_text(hjust = 0, size = 3) +
    xlim(1, 1.9) +
    labs(x = NULL, y = NULL, title = "P value") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 9))

  p_mean + p_diff + p_p + plot_layout(widths = c(1.2, 1.55, 0.42))
}
