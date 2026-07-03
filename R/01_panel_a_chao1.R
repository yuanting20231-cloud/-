source("R/00_config.R")
source("R/00_helpers.R")

make_panel_a <- function(dat = load_microbiome_data()) {
  richness <- as.data.frame(t(vegan::estimateR(dat$sample_taxa))) %>%
    rownames_to_column("Sample") %>%
    transmute(Sample, Chao1 = S.chao1) %>%
    left_join(dat$metadata, by = "Sample")

  plots <- lapply(levels(dat$metadata$Compartment), function(compartment_name) {
    d <- richness %>% filter(Compartment == compartment_name)
    letters <- dunn_letters(d, response = "Chao1", group = "Treatment")
    y_span <- diff(range(d$Chao1, na.rm = TRUE))
    y_pad <- ifelse(is.finite(y_span) && y_span > 0, 0.08 * y_span, 0.1)
    ann <- d %>%
      group_by(Treatment) %>%
      summarise(y = max(Chao1, na.rm = TRUE) + y_pad, .groups = "drop") %>%
      mutate(letter = letters[as.character(Treatment)])

    ggplot(d, aes(Treatment, Chao1, colour = Treatment)) +
      geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.55, fill = "white") +
      geom_jitter(width = 0.09, size = 2.3, alpha = 0.95) +
      geom_text(data = ann, aes(y = y, label = letter), inherit.aes = FALSE, size = 3.3) +
      scale_colour_manual(values = TREATMENT_COLORS, labels = TREATMENT_LABELS) +
      scale_x_discrete(labels = TREATMENT_LABELS) +
      labs(x = NULL, y = "Chao1 index", title = compartment_name) +
      theme_reference(10) +
      theme(legend.position = "none")
  })

  panel <- wrap_plots(plots, nrow = 1) + plot_annotation(tag_levels = NULL)
  save_plot(panel, "panel_a_chao1.png", FIG_WIDTHS$panel_a, 3.5)
  save_plot(panel, "panel_a_chao1.pdf", FIG_WIDTHS$panel_a, 3.5)
  panel
}

if (sys.nframe() == 0) make_panel_a()
