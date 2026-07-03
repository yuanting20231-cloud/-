source("R/00_config.R")
source("R/00_helpers.R")

make_panel_c <- function(dat = load_microbiome_data()) {
  family_long <- family_relative_long(dat)
  top <- top_families(family_long, n = TOP_N_FAMILIES)

  plot_data <- family_long %>%
    mutate(Family_plot = if_else(Family %in% top, Family, "Others")) %>%
    group_by(Group, Treatment, Compartment, Family_plot) %>%
    summarise(RelAbund = mean(RelAbund), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(RelAbund = 100 * RelAbund / sum(RelAbund)) %>%
    ungroup()

  family_order <- c("Others", rev(top))
  plot_data <- plot_data %>%
    mutate(
      Family_plot = factor(Family_plot, levels = family_order),
      Group = factor(Group, levels = dat$metadata %>% distinct(Group) %>% pull(Group))
    )

  palette <- c("grey70", grDevices::hcl.colors(length(top), palette = "Spectral", rev = TRUE))
  names(palette) <- family_order

  panel <- ggplot(plot_data, aes(Group, RelAbund, fill = Family_plot)) +
    geom_col(width = 0.72, colour = "grey35", linewidth = 0.16) +
    facet_grid(. ~ Compartment, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette, breaks = family_order, labels = family_order) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, 100)) +
    labs(x = NULL, y = "Percentage (%)", fill = "Family") +
    theme_reference(9) +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(0.35, "cm")
    ) +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE))

  save_plot(panel, "panel_c_family_composition.png", FIG_WIDTHS$panel_c, 6.2)
  save_plot(panel, "panel_c_family_composition.pdf", FIG_WIDTHS$panel_c, 6.2)
  panel
}

if (sys.nframe() == 0) make_panel_c()
