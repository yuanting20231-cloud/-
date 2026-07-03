source("R/00_config.R")
source("R/00_helpers.R")

make_panel_d <- function(dat = load_microbiome_data()) {
  family_long <- family_relative_long(dat)

  plots <- lapply(seq_len(nrow(COMPARISONS)), function(i) {
    g1 <- COMPARISONS$group1[i]
    g2 <- COMPARISONS$group2[i]
    prepared <- prepare_family_comparison(family_long, "Rhizosphere", g1, g2)
    plot_family_comparison(prepared) +
      plot_annotation(title = paste("Rhizosphere:", g1, "vs", g2)) &
      theme(plot.title = element_text(size = 10, hjust = 0))
  })

  panel <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots)))
  save_plot(panel, "panel_d_rhizosphere_differential.png", FIG_WIDTHS$panel_de, 7.0)
  save_plot(panel, "panel_d_rhizosphere_differential.pdf", FIG_WIDTHS$panel_de, 7.0)
  panel
}

if (sys.nframe() == 0) make_panel_d()
