# One-off: regenerate just the gam_monthly_dominance_<zone> GAM composite
# (plot_gam_figure()) after the plot_gam_figure()/gam_partial_effect() code
# changes, without
# redoing the slow GLM/regime/RF steps in run_full_analysis() (~57 min for
# those alone, unneeded here since none of that code changed).
# Run from repo root: Rscript func/render_gam_figure_only.R
source("func/driver_interactions.R")

render_gam_only <- function(plume_dir, fig_path){
  driver_matrices <- purrr::map(zones, ~ build_driver_matrix(.x, plume_dir = plume_dir)) |>
    purrr::set_names(zones)
  gam_models <- purrr::map(zones, ~ fit_gam(driver_matrices[[.x]])) |>
    purrr::set_names(zones) |> purrr::compact()
  plot_gam_figure(gam_models, fig_path)
  message("Wrote ", fig_path)
}

message("== GAM figure: dynamic threshold ==")
render_gam_only("output/panache/dynamic", "figures/ARTICLE/FIGURE_S7/Figure_S7.png")

message("== GAM figure: static threshold ==")
render_gam_only("output/panache/static", "figures/ARTICLE/FIGURE_S7/Figure_S7_static.png")
