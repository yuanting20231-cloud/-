# Figure 3a: change in network connectivity after random node removal.

make_fig3a_network_robustness <- function(dat) {
  counts <- family_count_matrix(dat)
  simulation_list <- list()
  network_list <- list()
  selected_taxa_list <- list()
  index <- 1L

  for (compartment in c("Rhizosphere", "Root")) {
    compartment_samples <- dat$metadata %>%
      filter(Compartment == compartment) %>%
      pull(Sample)

    taxa <- select_network_taxa(
      counts,
      sample_ids = compartment_samples,
      top_n = NETWORK_TOP_N,
      min_prevalence = NETWORK_MIN_PREVALENCE
    )

    selected_taxa_list[[compartment]] <- tibble(
      Compartment = compartment,
      Taxon = taxa,
      Rank = seq_along(taxa)
    )

    for (treatment in TREATMENT_LEVELS) {
      sample_ids <- dat$metadata %>%
        filter(Compartment == compartment, Treatment == treatment) %>%
        pull(Sample)

      graph <- build_association_network(
        counts[sample_ids, , drop = FALSE],
        taxa = taxa,
        group_label = paste(compartment, treatment)
      )

      simulation <- simulate_random_removal(
        graph,
        proportions = ROBUSTNESS_REMOVAL_PROPORTIONS,
        iterations = ROBUSTNESS_ITERATIONS,
        seed = RANDOM_SEED + index
      ) %>%
        mutate(
          Compartment = compartment,
          Treatment = treatment,
          Group = paste(compartment, treatment, sep = "__")
        )

      simulation_list[[index]] <- simulation
      network_list[[index]] <- tibble(
        Compartment = compartment,
        Treatment = treatment,
        Samples = length(sample_ids),
        Nodes = igraph::vcount(graph),
        Edges = igraph::ecount(graph),
        Density = igraph::edge_density(graph, loops = FALSE),
        MeanDegree = mean(igraph::degree(graph)),
        MeanClustering = suppressWarnings(igraph::transitivity(graph, type = "average"))
      )
      index <- index + 1L
    }
  }

  robustness <- bind_rows(simulation_list)
  network_summary <- bind_rows(network_list)
  selected_taxa <- bind_rows(selected_taxa_list)

  model_stats <- lapply(split(robustness, robustness$Group), function(d) {
    fit <- stats::lm(Connectivity ~ RemovalRatio, data = d)
    sm <- summary(fit)
    tibble(
      Compartment = d$Compartment[1],
      Treatment = d$Treatment[1],
      Intercept = unname(stats::coef(fit)[1]),
      Slope = unname(stats::coef(fit)[2]),
      AdjustedR2 = unname(sm$adj.r.squared),
      PValue = unname(sm$coefficients[2, 4])
    )
  }) %>% bind_rows() %>%
    mutate(
      Treatment = factor(Treatment, levels = TREATMENT_LEVELS),
      Label = paste0(
        "R² = ", formatC(AdjustedR2, format = "f", digits = 2),
        ", ", vapply(PValue, format_p_value, character(1))
      )
    )

  label_positions <- robustness %>%
    group_by(Compartment) %>%
    summarise(YMax = max(Connectivity, na.rm = TRUE), .groups = "drop")

  model_stats <- model_stats %>%
    left_join(label_positions, by = "Compartment") %>%
    group_by(Compartment) %>%
    arrange(Treatment, .by_group = TRUE) %>%
    mutate(
      LabelRank = row_number(),
      LabelY = ifelse(YMax > 0, YMax * (1.02 - 0.11 * (LabelRank - 1)), 0.01),
      LabelX = 0.98
    ) %>%
    ungroup()

  plot_data <- robustness %>%
    mutate(Treatment = factor(Treatment, levels = TREATMENT_LEVELS))

  plot <- ggplot(plot_data, aes(x = RemovalRatio, y = Connectivity, colour = Treatment)) +
    geom_point(
      position = position_jitter(width = 0.014, height = 0),
      alpha = 0.46,
      size = 1.25
    ) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7, alpha = 0.14) +
    geom_text(
      data = model_stats,
      aes(x = LabelX, y = LabelY, label = Label, colour = Treatment),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 2.55,
      show.legend = FALSE
    ) +
    facet_grid(rows = vars(Compartment), scales = "free_y") +
    scale_colour_manual(
      values = TREATMENT_COLORS,
      breaks = TREATMENT_LEVELS,
      labels = TREATMENT_LABELS
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      limits = c(0, 1),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +
    labs(
      x = "Ratio of randomly removed nodes",
      y = "Network connectivity",
      colour = NULL
    ) +
    theme_reference(9) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0, face = "plain", size = 9),
      panel.spacing.y = grid::unit(0.8, "lines"),
      plot.margin = margin(8, 8, 8, 8)
    )

  readr::write_csv(robustness, file.path(OUTPUT_DIR, "Figure3a_random_removal_raw.csv"))
  readr::write_csv(model_stats, file.path(OUTPUT_DIR, "Figure3a_linear_models.csv"))
  readr::write_csv(network_summary, file.path(OUTPUT_DIR, "Figure3a_network_summary.csv"))
  readr::write_csv(selected_taxa, file.path(OUTPUT_DIR, "Figure3a_selected_families.csv"))

  list(
    plot = plot,
    robustness = robustness,
    model_stats = model_stats,
    network_summary = network_summary,
    selected_taxa = selected_taxa
  )
}
