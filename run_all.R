source("R/00_config.R")
source("R/00_helpers.R")
source("R/01_panel_a_chao1.R")
source("R/02_panel_b_pcoa.R")
source("R/03_panel_c_family_composition.R")
source("R/04_panel_d_rhizosphere_differential.R")
source("R/05_panel_e_root.R")

dat <- load_microbiome_data()

panel_a <- make_panel_a(dat)
panel_b <- make_panel_b(dat)
panel_c <- make_panel_c(dat)
panel_d <- make_panel_d(dat)
panel_e <- make_panel_e(dat)

combined <- (panel_a | panel_b) /
  (panel_c | (panel_d / panel_e)) +
  plot_layout(heights = c(0.32, 0.68), widths = c(0.42, 0.58)) +
  plot_annotation(tag_levels = "a")

save_plot(combined, "Figure2_combined.png", FIG_WIDTHS$combined, 13.5)
save_plot(combined, "Figure2_combined.pdf", FIG_WIDTHS$combined, 13.5)

message("Finished. Results are in: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE))
