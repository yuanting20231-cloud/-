source("R/00_config.R")
source("R/00_helpers.R")
source("R/06_network_helpers.R")
source("R/07_fig3a_network_robustness.R")
source("R/08_fig3b_netshift.R")
source("R/09_fig3c_driver_upset.R")

if (!file.exists(PATHS$taxonomy)) {
  stop(
    "Figure 3 requires data/taxonomy.tsv with columns '#ID' and 'family'.",
    call. = FALSE
  )
}

dat <- load_microbiome_data()

figure3a <- make_fig3a_network_robustness(dat)
figure3b <- make_fig3b_netshift(dat)
figure3c <- make_fig3c_driver_upset(figure3b)

save_plot(
  figure3a$plot,
  "Figure3a_network_robustness.png",
  FIG3_WIDTHS$panel_a,
  7.4
)
save_plot(
  figure3a$plot,
  "Figure3a_network_robustness.pdf",
  FIG3_WIDTHS$panel_a,
  7.4
)

save_plot(
  figure3b$plot,
  "Figure3b_netshift.png",
  FIG3_WIDTHS$panel_b,
  9.0
)
save_plot(
  figure3b$plot,
  "Figure3b_netshift.pdf",
  FIG3_WIDTHS$panel_b,
  9.0
)

save_plot(
  figure3c$plot,
  "Figure3c_driver_upset.png",
  FIG3_WIDTHS$panel_c,
  4.8
)
save_plot(
  figure3c$plot,
  "Figure3c_driver_upset.pdf",
  FIG3_WIDTHS$panel_c,
  4.8
)

left_column <- (figure3a$plot / figure3c$plot) +
  patchwork::plot_layout(heights = c(0.70, 0.30))
combined <- (left_column | figure3b$plot) +
  patchwork::plot_layout(widths = c(0.34, 0.66))

save_plot(
  combined,
  "Figure3_combined.png",
  FIG3_WIDTHS$combined,
  10.8
)
save_plot(
  combined,
  "Figure3_combined.pdf",
  FIG3_WIDTHS$combined,
  10.8
)

capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "Figure3_sessionInfo.txt"))
message(
  "Figure 3 finished. Results are in: ",
  normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)
)
