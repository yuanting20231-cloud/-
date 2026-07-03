source("R/00_config.R")
source("R/00_helpers.R")

make_pcoa_one <- function(dat, compartment_name) {
  meta <- dat$metadata %>% filter(Compartment == compartment_name) %>% droplevels()
  rel <- relative_abundance_matrix(dat$sample_taxa[meta$Sample, , drop = FALSE])
  bray <- vegan::vegdist(rel, method = "bray")
  ord <- vegan::wcmdscale(bray, k = 2, eig = TRUE, add = "lingoes")

  coords <- as.data.frame(ord$points[, 1:2, drop = FALSE])
  names(coords) <- c("PCoA1", "PCoA2")
  coords <- coords %>%
    rownames_to_column("Sample") %>%
    left_join(meta, by = "Sample")

  positive_eig <- ord$eig[ord$eig > 0]
  if (length(positive_eig) < 2) stop("PCoA returned fewer than two positive eigenvalues.")
  explained <- 100 * positive_eig[1:2] / sum(positive_eig)

  set.seed(RANDOM_SEED)
  perm <- vegan::adonis2(bray ~ Treatment, data = meta, permutations = PERMUTATIONS)
  r2 <- perm$R2[1]
  p_value <- perm$`Pr(>F)`[1]
  p_text <- if (is.na(p_value)) "NA" else if (p_value < 0.001) "< 0.001" else paste0("= ", format(round(p_value, 3), nsmall = 3))
  stat_label <- sprintf("R² = %.2f  P %s", r2, p_text)

  ggplot(coords, aes(PCoA1, PCoA2, colour = Treatment, fill = Treatment)) +
    stat_ellipse(type = "t", level = 0.95, linewidth = 0.9, alpha = 0, show.legend = FALSE) +
    geom_point(shape = 21, size = 3, stroke = 0.55) +
    annotate("text", x = Inf, y = Inf, label = stat_label,
             hjust = 1.05, vjust = 1.2, size = 3.2) +
    scale_colour_manual(values = TREATMENT_COLORS, labels = TREATMENT_LABELS) +
    scale_fill_manual(values = TREATMENT_COLORS, labels = TREATMENT_LABELS) +
    labs(
      x = sprintf("PCoA 1 (%.2f%%)", explained[1]),
      y = sprintf("PCoA 2 (%.2f%%)", explained[2]),
      title = compartment_name
    ) +
    theme_reference(10) +
    theme(legend.position = c(0.86, 0.18))
}

make_panel_b <- function(dat = load_microbiome_data()) {
  plots <- lapply(levels(dat$metadata$Compartment), function(x) make_pcoa_one(dat, x))
  panel <- wrap_plots(plots, nrow = 1)
  save_plot(panel, "panel_b_pcoa.png", FIG_WIDTHS$panel_b, 3.5)
  save_plot(panel, "panel_b_pcoa.pdf", FIG_WIDTHS$panel_b, 3.5)
  panel
}

if (sys.nframe() == 0) make_panel_b()
