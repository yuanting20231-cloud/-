# Figure 3c: UpSet-style intersection plot for driver taxa.

make_fig3c_driver_upset <- function(netshift_result) {
  drivers <- netshift_result$drivers %>%
    transmute(
      Taxon,
      Comparison,
      Compartment,
      Set = paste0(Comparison, " (", Compartment, ")")
    ) %>%
    distinct()

  comparison_order <- NETSHIFT_COMPARISONS$comparison_id
  set_order <- c(
    paste0(comparison_order, " (Root)"),
    paste0(comparison_order, " (Rhizosphere)")
  )
  set_order <- set_order[set_order %in% unique(drivers$Set)]

  if (nrow(drivers) == 0 || length(set_order) == 0) {
    stop("No driver taxa were available for the UpSet plot.", call. = FALSE)
  }

  all_taxa <- sort(unique(drivers$Taxon))
  membership_long <- tidyr::crossing(Taxon = all_taxa, Set = set_order) %>%
    left_join(drivers %>% mutate(Present = TRUE), by = c("Taxon", "Set")) %>%
    mutate(Present = dplyr::coalesce(Present, FALSE))

  membership_wide <- membership_long %>%
    select(Taxon, Set, Present) %>%
    pivot_wider(names_from = Set, values_from = Present, values_fill = FALSE)

  binary_matrix <- as.matrix(membership_wide[, set_order, drop = FALSE])
  storage.mode(binary_matrix) <- "integer"
  taxon_patterns <- tibble(
    Taxon = membership_wide$Taxon,
    Pattern = apply(binary_matrix, 1, paste0, collapse = "")
  ) %>%
    filter(grepl("1", Pattern, fixed = TRUE))

  combination_table_all <- taxon_patterns %>%
    group_by(Pattern) %>%
    summarise(
      IntersectionSize = n(),
      Taxa = paste(sort(Taxon), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(NumberOfSets = stringr::str_count(Pattern, "1")) %>%
    arrange(desc(IntersectionSize), desc(NumberOfSets), Pattern)

  keep_n <- min(UPSET_MAX_COMBINATIONS, nrow(combination_table_all))
  combination_table <- combination_table_all %>%
    slice_head(n = keep_n) %>%
    mutate(
      ComboIndex = row_number(),
      Combination = paste0("C", ComboIndex)
    )

  if (nrow(combination_table) == 0) {
    stop("No non-empty driver-taxa intersections were found.", call. = FALSE)
  }

  matrix_rows <- lapply(seq_len(nrow(combination_table)), function(i) {
    bits <- strsplit(combination_table$Pattern[i], "", fixed = TRUE)[[1]]
    tibble(
      ComboIndex = combination_table$ComboIndex[i],
      Combination = combination_table$Combination[i],
      Set = set_order,
      Present = bits == "1"
    )
  }) %>% bind_rows()

  display_set_order <- rev(set_order)
  matrix_rows <- matrix_rows %>%
    mutate(SetY = match(Set, display_set_order))

  connection_segments <- matrix_rows %>%
    filter(Present) %>%
    group_by(ComboIndex) %>%
    summarise(
      YMin = min(SetY),
      YMax = max(SetY),
      .groups = "drop"
    )

  set_sizes <- membership_long %>%
    filter(Present) %>%
    count(Set, name = "SetSize") %>%
    right_join(tibble(Set = set_order), by = "Set") %>%
    mutate(
      SetSize = dplyr::coalesce(SetSize, 0L),
      SetY = match(Set, display_set_order)
    )

  multi_set <- combination_table %>% filter(NumberOfSets > 1)
  highlight_index <- if (nrow(multi_set) > 0) {
    multi_set$ComboIndex[which.max(multi_set$IntersectionSize)]
  } else {
    combination_table$ComboIndex[which.max(combination_table$IntersectionSize)]
  }

  combination_table <- combination_table %>%
    mutate(Highlight = ComboIndex == highlight_index)

  highlighted <- combination_table %>% filter(Highlight) %>% slice_head(n = 1)
  highlighted_taxa <- strsplit(highlighted$Taxa, "; ", fixed = TRUE)[[1]]
  highlighted_label <- paste(head(highlighted_taxa, 5), collapse = "\n")
  if (length(highlighted_taxa) > 5) {
    highlighted_label <- paste0(highlighted_label, "\n...")
  }

  max_intersection <- max(combination_table$IntersectionSize)
  p_top <- ggplot(combination_table, aes(x = ComboIndex, y = IntersectionSize)) +
    geom_col(aes(fill = Highlight), width = 0.68) +
    scale_fill_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = "#D73027"), guide = "none") +
    geom_text(aes(label = IntersectionSize), vjust = -0.35, size = 3) +
    geom_label(
      data = highlighted,
      aes(
        x = ComboIndex,
        y = IntersectionSize + max(0.6, max_intersection * 0.18)
      ),
      label = highlighted_label,
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0,
      size = 2.6,
      label.size = 0.2,
      fill = "white"
    ) +
    scale_x_continuous(
      limits = c(0.45, max(combination_table$ComboIndex) + 0.55),
      breaks = combination_table$ComboIndex,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max_intersection * 1.55 + 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, y = "Driver taxa in each set") +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(6, 8, 0, 0)
    )

  p_matrix <- ggplot(matrix_rows, aes(x = ComboIndex, y = SetY)) +
    geom_segment(
      data = connection_segments,
      aes(x = ComboIndex, xend = ComboIndex, y = YMin, yend = YMax),
      inherit.aes = FALSE,
      linewidth = 0.55,
      colour = "black"
    ) +
    geom_point(aes(fill = Present), shape = 21, size = 2.5, stroke = 0.25) +
    scale_fill_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "black"), guide = "none") +
    scale_x_continuous(
      limits = c(0.45, max(combination_table$ComboIndex) + 0.55),
      breaks = combination_table$ComboIndex,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(display_set_order),
      labels = display_set_order,
      limits = c(0.5, length(display_set_order) + 0.5),
      position = "right"
    ) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 8.5) +
    theme(
      axis.text.y = element_text(colour = "black", hjust = 0, margin = margin(l = 5)),
      axis.text.x = element_blank(),
      plot.margin = margin(0, 75, 6, 0)
    )

  p_set <- ggplot(set_sizes, aes(x = SetSize, y = SetY)) +
    geom_col(
      orientation = "y",
      width = 0.62,
      fill = "#BDBDBD",
      colour = "#7F7F7F",
      linewidth = 0.25
    ) +
    geom_text(aes(label = SetSize), hjust = 1.15, size = 2.7) +
    scale_x_reverse(expand = expansion(mult = c(0.08, 0.18))) +
    scale_y_continuous(
      breaks = seq_along(display_set_order),
      labels = NULL,
      limits = c(0.5, length(display_set_order) + 0.5)
    ) +
    labs(x = "Set Size", y = NULL) +
    theme_classic(base_size = 8.5) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 4, 6, 4)
    )

  plot <- (patchwork::plot_spacer() | p_top) /
    (p_set | p_matrix) +
    patchwork::plot_layout(widths = c(0.28, 1), heights = c(0.68, 0.32))

  readr::write_csv(
    membership_long,
    file.path(OUTPUT_DIR, "Figure3c_driver_set_membership.csv")
  )
  readr::write_csv(
    combination_table,
    file.path(OUTPUT_DIR, "Figure3c_driver_intersections.csv")
  )

  list(
    plot = plot,
    membership = membership_long,
    intersections = combination_table,
    set_sizes = set_sizes
  )
}
