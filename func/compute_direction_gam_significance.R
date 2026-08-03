# One-off: extract the joint (all-octant) significance of the categorical
# wind_dir_cat / wave_dir_cat parametric terms from the multi-driver GAM, for
# manuscript Table 6's "Wind direction" / "Wave direction" rows. fit_gam()'s
# tensor-product smooths only take numeric driver pairs, so wind_dir_cat and
# wave_dir_cat enter as flat parametric main-effect terms instead
# (func/driver_interactions.R::fit_gam()) -- summary.gam()'s $s.table (what
# output/STATS/driver_gam_summary.csv already saves) only covers the te()
# smooth terms, not these parametric ones. summary.gam()'s $pTerms.table is
# the right extraction: unlike $p.table (one row per dummy-coded factor
# level, i.e. per compass octant), $pTerms.table reports one joint
# significance test per whole parametric term (e.g. one row for
# "wind_dir_cat" testing all 8 octants together).
#
# Reuses build_driver_matrix()/fit_gam() directly rather than
# Apply_driver_interactions_analysis()'s full run_full_analysis() pipeline,
# since only the (already-fitted, already fast) GAM step is needed here --
# not the slow random-forest H-statistic step that dominated the ~57 min
# full rerun.
#
# Run from repo root: Rscript func/compute_direction_gam_significance.R
source("func/driver_interactions.R")

extract_direction_terms <- function(plume_dir, stats_dir){
  driver_matrices <- purrr::map(zones, ~ build_driver_matrix(.x, plume_dir = plume_dir)) |>
    purrr::set_names(zones)

  gam_models <- purrr::map(zones, ~ fit_gam(driver_matrices[[.x]])) |>
    purrr::set_names(zones) |> purrr::compact()

  direction_terms <- purrr::imap_dfr(gam_models, function(m, zone_name){
    s <- summary(m)
    pt <- s$pTerms.table
    if(is.null(pt)) return(NULL)
    pt <- pt[rownames(pt) %in% c("wind_dir_cat", "wave_dir_cat"), , drop = FALSE]
    if(nrow(pt) == 0) return(NULL)
    tibble::tibble(zone = zone_name, term = rownames(pt),
                   df = pt[, "df"], F = pt[, "F"], p = pt[, "p-value"])
  })

  readr::write_csv(direction_terms, file.path(stats_dir, "driver_gam_direction_terms.csv"))
  direction_terms
}

message("== Direction GAM significance: dynamic threshold ==")
dynamic_terms <- extract_direction_terms("output/panache/dynamic", "output/STATS")
print(dynamic_terms, n = Inf)

message("== Direction GAM significance: static threshold ==")
static_terms <- extract_direction_terms("output/panache/static", "output/STATS/static")
print(static_terms, n = Inf)
