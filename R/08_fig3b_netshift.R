# Figure 3b: NetShift-style circular comparison networks.

make_fig3b_netshift <- function(dat) {
  counts <- family_count_matrix(dat)
  results <- list()
  plot_list <- list()
  selected_taxa_list <- list()
  network_summary_list <- list()
  index <- 1L

  for (compartment in c("Rhizosphere", "Root")) {
    for (i in seq_len(nrow(NETSHIFT_COMPARISONS))) {
      comparison <- NETSHIFT_COMPARISONS[i, ]
      case_samples <- dat$metadata %>%
        filter(Compartment == compartment, Treatment == comparison$case) %>%
        pull(Sample)
      control_samples <- dat$metadata %>%
        filter(Compartment == compartment, Treatment == comparison$control) %>%
        pull(Sample)

      taxa <- select_network_taxa(
        counts,
        sample_ids = c(case_samples, control_samples),
        top_n = NETWORK_TOP_N,
        min_prevalence = NETWORK_MIN_PREVALENCE,
        require_present_in = list(case_samples, control_samples)
      )

      selected_taxa_list[[index]] <- tibble(
        Comparison = comparison$comparison_id,
        Compartment = compartment,
        Taxon = taxa,
        Rank = seq_along(taxa)
      )

      control_graph <- build_association_network(
        counts[control_samples, , drop = FALSE],
        taxa = taxa,
        group_label = paste(compartment, comparison$control_label, "control")
      )
      case_graph <- build_association_network(
        counts[case_samples, , drop = FALSE],
        taxa = taxa,
        group_label = paste(compartment, comparison$case_label, "case")
      )

      result <- calculate_netshift(
        control_graph = control_graph,
        case_graph = case_graph,
        comparison_id = comparison$comparison_id,
        compartment = compartment,
        control_label = comparison$control_label,
        case_label = comparison$case_label
      )

      edge_union_n <- nrow(result$edges)
      edge_both_n <- sum(result$edges$edge_class == "Both")
      jei <- if (edge_union_n == 0) NA_real_ else edge_both_n / edge_union_n

      network_summary_list[[index]] <- tibble(
        Comparison = comparison$comparison_id,
        Compartment = compartment,
        Control = comparison$control_label,
        Case = comparison$case_label,
        CommonNodes = igraph::vcount(result$control_graph),
        ControlEdges = igraph::ecount(result$control_graph),
        CaseEdges = igraph::ecount(result$case_graph),
        SharedEdges = edge_both_n,
        UnionEdges = edge_union_n,
        JaccardEdgeIndex = jei,
        DriverTaxa = sum(result$nodes$Driver)
      )

      results[[index]] <- result
      plot_list[[index]] <- plot_netshift_network(result, show_legend = TRUE)
      index <- index + 1L
    }
  }

  combined_plot <- patchwork::wrap_plots(plot_list, ncol = nrow(NETSHIFT_COMPARISONS), guides = "collect") &
    theme(legend.position = "bottom")

  node_table <- bind_rows(lapply(results, function(x) x$nodes))
  edge_table <- bind_rows(lapply(results, function(x) {
    x$edges %>% mutate(Comparison = x$comparison_id, Compartment = x$compartment)
  }))
  selected_taxa <- bind_rows(selected_taxa_list)
  network_summary <- bind_rows(network_summary_list)
  driver_table <- node_table %>%
    filter(Driver) %>%
    arrange(Compartment, Comparison, desc(NESH), desc(DeltaBetweenness))

  readr::write_csv(node_table, file.path(OUTPUT_DIR, "Figure3b_netshift_node_statistics.csv"))
  readr::write_csv(edge_table, file.path(OUTPUT_DIR, "Figure3b_netshift_edge_classes.csv"))
  readr::write_csv(driver_table, file.path(OUTPUT_DIR, "Figure3b_driver_taxa.csv"))
  readr::write_csv(network_summary, file.path(OUTPUT_DIR, "Figure3b_network_comparison_summary.csv"))
  readr::write_csv(selected_taxa, file.path(OUTPUT_DIR, "Figure3b_selected_families.csv"))

  list(
    plot = combined_plot,
    results = results,
    nodes = node_table,
    edges = edge_table,
    drivers = driver_table,
    network_summary = network_summary,
    selected_taxa = selected_taxa
  )
}
