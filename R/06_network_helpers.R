# Shared functions for Figure 3 network robustness and NetShift-style analyses.

required_network_packages <- c("igraph", "tidygraph", "ggraph", "corpcor", "scales")
missing_network_packages <- required_network_packages[
  !vapply(required_network_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_network_packages) > 0) {
  stop(
    "Missing network packages: ", paste(missing_network_packages, collapse = ", "),
    "\nRun source('install_packages.R') first.",
    call. = FALSE
  )
}

safe_scale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

family_count_matrix <- function(dat) {
  if (!"family" %in% names(dat$taxonomy)) {
    stop("taxonomy.tsv must contain a 'family' column for Figure 3.", call. = FALSE)
  }

  tax_family <- dat$taxonomy %>%
    transmute(
      `#ID`,
      Family = clean_taxon(family, fallback = "Unclassified")
    )

  family_long <- dat$abundance %>%
    pivot_longer(-`#ID`, names_to = "Sample", values_to = "Count") %>%
    left_join(tax_family, by = "#ID") %>%
    mutate(Family = clean_taxon(Family, fallback = "Unclassified")) %>%
    group_by(Sample, Family) %>%
    summarise(Count = sum(Count), .groups = "drop")

  family_wide <- family_long %>%
    tidyr::pivot_wider(names_from = Family, values_from = Count, values_fill = 0) %>%
    arrange(match(Sample, dat$metadata$Sample))

  mat <- family_wide %>%
    tibble::column_to_rownames("Sample") %>%
    as.matrix()
  storage.mode(mat) <- "numeric"
  mat[dat$metadata$Sample, , drop = FALSE]
}

select_network_taxa <- function(count_matrix, sample_ids, top_n = NETWORK_TOP_N,
                                min_prevalence = NETWORK_MIN_PREVALENCE,
                                require_present_in = NULL) {
  x <- count_matrix[sample_ids, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 0) stop("No samples or taxa available for network selection.", call. = FALSE)

  prevalence <- colMeans(x > 0)
  rel <- sweep(x, 1, pmax(rowSums(x), 1), "/")
  mean_rel <- colMeans(rel)
  keep <- prevalence >= min_prevalence & colSums(x) > 0

  if (!is.null(require_present_in)) {
    for (ids in require_present_in) {
      keep <- keep & colSums(count_matrix[ids, , drop = FALSE] > 0) > 0
    }
  }

  ranked <- tibble(
    Taxon = colnames(x),
    Prevalence = prevalence,
    MeanRelativeAbundance = mean_rel,
    Keep = keep
  ) %>%
    filter(Keep) %>%
    arrange(desc(MeanRelativeAbundance), desc(Prevalence), Taxon)

  if (nrow(ranked) < 3) {
    stop("Fewer than three taxa passed the network filtering criteria.", call. = FALSE)
  }

  head(ranked$Taxon, min(top_n, nrow(ranked)))
}

clr_transform <- function(count_matrix, pseudocount = NETWORK_PSEUDOCOUNT) {
  log_x <- log(count_matrix + pseudocount)
  sweep(log_x, 1, rowMeans(log_x), "-")
}

build_association_network <- function(count_matrix, taxa, group_label = "network") {
  x <- count_matrix[, taxa, drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]

  if (nrow(x) < 6) {
    warning(
      group_label, " has only ", nrow(x),
      " samples. Shrinkage correlations will run, but network topology is exploratory."
    )
  }
  if (ncol(x) < 3) stop(group_label, " has fewer than three taxa.", call. = FALSE)

  clr <- clr_transform(x)
  variable_taxa <- apply(clr, 2, stats::sd, na.rm = TRUE) > 0
  clr <- clr[, variable_taxa, drop = FALSE]
  if (ncol(clr) < 3) stop(group_label, " has fewer than three variable taxa.", call. = FALSE)

  cor_mat <- suppressWarnings(corpcor::cor.shrink(clr, verbose = FALSE))
  diag(cor_mat) <- 0

  idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  edge_candidates <- tibble(
    from = colnames(cor_mat)[idx[, 1]],
    to = colnames(cor_mat)[idx[, 2]],
    correlation = cor_mat[idx],
    abs_correlation = abs(cor_mat[idx]),
    association_sign = ifelse(cor_mat[idx] >= 0, "Positive", "Negative")
  ) %>%
    filter(is.finite(correlation))

  positive_candidates <- edge_candidates %>% filter(correlation > 0)
  candidate_pool <- if (isTRUE(NETWORK_POSITIVE_ONLY) && nrow(positive_candidates) > 0) {
    positive_candidates
  } else {
    if (isTRUE(NETWORK_POSITIVE_ONLY)) {
      warning(group_label, " had no positive correlations; selecting by absolute correlation.")
    }
    edge_candidates
  }

  candidate_pool <- candidate_pool %>% arrange(desc(abs_correlation))
  if (nrow(candidate_pool) == 0) {
    graph <- igraph::make_empty_graph(n = ncol(clr), directed = FALSE)
    igraph::V(graph)$name <- colnames(clr)
    return(graph)
  }

  q_cut <- as.numeric(stats::quantile(
    candidate_pool$abs_correlation,
    probs = NETWORK_EDGE_QUANTILE,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  ))
  threshold <- max(NETWORK_MIN_ABS_COR, q_cut)
  selected <- candidate_pool %>% filter(abs_correlation >= threshold)

  if (nrow(selected) < NETWORK_MIN_EDGES) {
    selected <- candidate_pool %>% slice_head(n = min(NETWORK_MIN_EDGES, nrow(candidate_pool)))
  }
  if (nrow(selected) > NETWORK_MAX_EDGES) {
    selected <- selected %>% slice_head(n = NETWORK_MAX_EDGES)
  }

  vertices <- tibble(name = colnames(clr))
  graph <- igraph::graph_from_data_frame(
    selected %>% select(from, to, correlation, abs_correlation, association_sign),
    directed = FALSE,
    vertices = vertices
  )
  igraph::simplify(
    graph,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(
      correlation = "mean",
      abs_correlation = "max",
      association_sign = "first"
    )
  )
}

edge_key_table <- function(graph) {
  if (igraph::ecount(graph) == 0) {
    return(tibble(from = character(), to = character(), key = character()))
  }
  edges <- igraph::as_data_frame(graph, what = "edges") %>%
    transmute(
      from = pmin(as.character(from), as.character(to)),
      to = pmax(as.character(from), as.character(to))
    ) %>%
    distinct(from, to) %>%
    mutate(key = paste(from, to, sep = "||"))
  edges
}

induce_common_nodes <- function(graph_a, graph_b) {
  common <- intersect(igraph::V(graph_a)$name, igraph::V(graph_b)$name)
  if (length(common) < 3) stop("The two networks have fewer than three common taxa.", call. = FALSE)
  list(
    a = igraph::induced_subgraph(graph_a, vids = common),
    b = igraph::induced_subgraph(graph_b, vids = common)
  )
}

calculate_netshift <- function(control_graph, case_graph, comparison_id, compartment,
                               control_label, case_label) {
  common_graphs <- induce_common_nodes(control_graph, case_graph)
  control_graph <- common_graphs$a
  case_graph <- common_graphs$b
  nodes <- sort(intersect(igraph::V(control_graph)$name, igraph::V(case_graph)$name))

  control_edges <- edge_key_table(control_graph)
  case_edges <- edge_key_table(case_graph)

  union_edges <- full_join(
    control_edges %>% mutate(in_control = TRUE),
    case_edges %>% mutate(in_case = TRUE),
    by = c("from", "to", "key")
  ) %>%
    mutate(
      in_control = dplyr::coalesce(in_control, FALSE),
      in_case = dplyr::coalesce(in_case, FALSE),
      edge_class = case_when(
        in_control & in_case ~ "Both",
        in_control ~ "Control only",
        TRUE ~ "Case only"
      ),
      edge_class = factor(edge_class, levels = c("Control only", "Case only", "Both"))
    )

  degree_case <- igraph::degree(case_graph, v = nodes)
  max_degree_case <- max(c(degree_case, 1))

  control_bet <- igraph::betweenness(control_graph, v = nodes, directed = FALSE, normalized = TRUE)
  case_bet <- igraph::betweenness(case_graph, v = nodes, directed = FALSE, normalized = TRUE)
  control_bet_scaled <- safe_scale01(control_bet)
  case_bet_scaled <- safe_scale01(case_bet)

  node_stats <- lapply(seq_along(nodes), function(i) {
    node <- nodes[i]
    n_control <- as.character(igraph::neighbors(control_graph, node))
    n_case <- as.character(igraph::neighbors(case_graph, node))
    n_union <- union(n_control, n_case)
    n_intersection <- intersect(n_control, n_case)
    n_case_only <- setdiff(n_case, n_control)

    if (length(n_union) == 0) {
      x <- 1
      y <- 0
      z <- 0
    } else {
      x <- length(n_intersection) / length(n_union)
      y <- length(n_case_only) / max_degree_case
      z <- length(n_case_only) / length(n_union)
    }

    tibble(
      Taxon = node,
      DegreeControl = length(n_control),
      DegreeCase = length(n_case),
      SharedNeighbors = length(n_intersection),
      CaseOnlyNeighbors = length(n_case_only),
      UnionNeighbors = length(n_union),
      NESH = 1 - x + y + z,
      BetweennessControl = as.numeric(control_bet[i]),
      BetweennessCase = as.numeric(case_bet[i]),
      ScaledBetweennessControl = control_bet_scaled[i],
      ScaledBetweennessCase = case_bet_scaled[i],
      DeltaBetweenness = case_bet_scaled[i] - control_bet_scaled[i]
    )
  }) %>% bind_rows()

  node_stats <- node_stats %>%
    mutate(Driver = DeltaBetweenness > 0 & NESH >= NESH_DRIVER_MIN)

  if (!any(node_stats$Driver)) {
    fallback <- node_stats %>%
      filter(DeltaBetweenness > 0) %>%
      arrange(desc(NESH), desc(DeltaBetweenness)) %>%
      slice_head(n = min(NESH_FALLBACK_TOP_N, n())) %>%
      pull(Taxon)
    if (length(fallback) == 0) {
      fallback <- node_stats %>%
        arrange(desc(NESH), desc(DeltaBetweenness)) %>%
        slice_head(n = min(NESH_FALLBACK_TOP_N, n())) %>%
        pull(Taxon)
    }
    node_stats$Driver <- node_stats$Taxon %in% fallback
    warning(comparison_id, " ", compartment, ": no taxa met the formal driver rule; top taxa were highlighted as fallback.")
  }

  membership <- rep(1L, length(nodes))
  names(membership) <- nodes
  if (igraph::ecount(case_graph) > 0) {
    case_for_cluster <- case_graph
    igraph::E(case_for_cluster)$weight <- 1
    membership <- igraph::membership(igraph::cluster_louvain(case_for_cluster, weights = NULL))
  }

  node_stats <- node_stats %>%
    mutate(
      Community = as.integer(membership[Taxon]),
      Comparison = comparison_id,
      Compartment = compartment,
      Control = control_label,
      Case = case_label
    ) %>%
    arrange(Community, desc(Driver), desc(NESH), Taxon)

  vertex_table <- node_stats %>%
    transmute(
      name = Taxon,
      NESH,
      DeltaBetweenness,
      Driver,
      Community
    )

  union_graph <- igraph::graph_from_data_frame(
    union_edges %>% select(from, to, edge_class),
    directed = FALSE,
    vertices = vertex_table
  )

  list(
    comparison_id = comparison_id,
    compartment = compartment,
    control_label = control_label,
    case_label = case_label,
    control_graph = control_graph,
    case_graph = case_graph,
    union_graph = union_graph,
    nodes = node_stats,
    edges = union_edges
  )
}

simulate_random_removal <- function(graph, proportions = ROBUSTNESS_REMOVAL_PROPORTIONS,
                                    iterations = ROBUSTNESS_ITERATIONS,
                                    seed = RANDOM_SEED) {
  set.seed(seed)
  n0 <- igraph::vcount(graph)
  if (n0 < 3) stop("Network robustness requires at least three nodes.", call. = FALSE)
  denominator <- n0 * (n0 - 1)

  out <- lapply(proportions, function(prop) {
    remove_n <- min(floor(prop * n0), n0 - 1L)
    reps <- if (remove_n == 0) 1L else iterations

    lapply(seq_len(reps), function(iteration) {
      removed <- if (remove_n == 0) character() else sample(igraph::V(graph)$name, remove_n)
      reduced <- igraph::delete_vertices(graph, removed)
      components <- if (igraph::vcount(reduced) > 0) igraph::components(reduced)$csize else 0
      largest_component <- if (length(components) == 0) 0 else max(components)

      tibble(
        RemovalRatio = prop,
        Iteration = iteration,
        NodesRemaining = igraph::vcount(reduced),
        EdgesRemaining = igraph::ecount(reduced),
        Connectivity = if (denominator == 0) 0 else 2 * igraph::ecount(reduced) / denominator,
        LargestComponentFraction = largest_component / n0
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  out
}

format_p_value <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", formatC(p, format = "f", digits = 3))
}

plot_netshift_network <- function(result, show_legend = TRUE) {
  graph <- result$union_graph
  if (igraph::vcount(graph) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(result$compartment, result$case_label, "vs", result$control_label)))
  }

  layout <- ggraph::create_layout(graph, layout = "linear", circular = TRUE)
  n <- nrow(layout)
  angle <- 90 - 360 * (seq_len(n) - 1) / n
  hjust <- ifelse(angle < -90, 1, 0)
  angle <- ifelse(angle < -90, angle + 180, angle)
  layout$label_angle <- angle
  layout$label_hjust <- hjust
  layout$label_x <- layout$x * 1.16
  layout$label_y <- layout$y * 1.16

  p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_arc(
      aes(edge_colour = edge_class),
      strength = 0.86,
      edge_width = 0.35,
      alpha = 0.62,
      show.legend = show_legend
    ) +
    ggraph::scale_edge_colour_manual(
      values = NETSHIFT_EDGE_COLORS,
      drop = FALSE,
      name = NULL
    ) +
    ggraph::geom_node_point(
      aes(size = NESH, colour = Driver),
      alpha = 0.95,
      stroke = 0.25
    ) +
    scale_size_continuous(range = c(1.5, 5.5), name = "NESH") +
    scale_colour_manual(values = c(`TRUE` = "#D73027", `FALSE` = "#222222"), guide = "none") +
    geom_text(
      data = layout,
      aes(
        x = label_x,
        y = label_y,
        label = name,
        angle = label_angle,
        hjust = label_hjust,
        colour = Driver
      ),
      inherit.aes = FALSE,
      size = 2.25,
      lineheight = 0.9
    ) +
    annotate(
      "label",
      x = 0,
      y = 0,
      label = paste0(result$case_label, " vs ", result$control_label),
      size = 3,
      label.size = 0.2,
      fill = "white"
    ) +
    coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45), clip = "off") +
    labs(title = result$compartment) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
      plot.margin = margin(14, 26, 14, 26),
      legend.position = if (show_legend) "bottom" else "none"
    )

  p
}
