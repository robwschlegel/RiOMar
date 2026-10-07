# func/figure.R
# This is called by func/figure.py to provide graphical outputs


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(scales)
library(maps)
library(ggpubr)
library(cowplot)
library(magick)

# For zones_bbox
source("func/util.R")

# For zone_meta/get_zone_meta/combine_plume_driver/plot_driver_rose/etc
source("func/multi.R")

# For the tide QC diagnostics (tide_qc_all(), save_tide_qc_examples(), etc.)
# -- .load_tide_raw() and the QC logic load_tide_gauge() calls both live in
# func/util.R itself, not here; nothing in the pipeline path requires this
# file, it is sourced only for its standalone diagnostic/reporting functions.
source("func/tide.R")


# Utils -------------------------------------------------------------------

# High-resolution France coastline (GADM level-0 boundary, ~216,000 vertices)
# for the study-zone-map figure's zoomed insets.
# Loaded lazily (library(sf) is NOT called at file scope

# Step size for four_ticks_from_zero() below, factored out so a dual-axis
# figure (e.g. plot_plume_area_timeseries()) can also compute the step needed to
# reach an arbitrary target -- not just a series' own max(x) -- letting two
# independently-scaled axes share the same top tick instead of each rounding
# up separately and leaving whichever axis rounds up less stranded well
# below the panel top (fixed 2026-08-11, the plume-area-timeseries figure's
# "dead space" tweak).
# Unlike base pretty()/scales::breaks_pretty(), which restrict the step to a
# 1/2/5 x 10^k multiple, a step of 2.5 x 10^k is allowed too -- needed so 3
# equal steps can land close to the target (e.g. target ~7500 -> steps of
# 2500, giving 0, 2500, 5000, 7500) rather than overshooting it by a lot.
# Picks the *smallest* nice step whose 3 steps still reach the target --
# picking whichever nice step is merely closest to target/3 (tried first)
# can round down and clip real data above the top tick, found via the Gulf
# of Lion's actual max area (8987 km^2) rounding down to a 7500 top tick.
nice_step_for_target <- function(target){
  raw_step <- target / 3
  magnitude <- 10 ^ floor(log10(raw_step))
  candidate_steps <- c(1, 2, 2.5, 5, 10) * magnitude
  min(candidate_steps[candidate_steps * 3 >= target])
}

# Four evenly spaced y-axis ticks from 0 up to a "nice" round number at or
# above max(x). See nice_step_for_target() above for the step-selection logic.
four_ticks_from_zero <- function(x){
  step <- nice_step_for_target(max(x, na.rm = TRUE))
  seq(0, step * 3, by = step)
}

# Three evenly spaced y-axis ticks that avoid landing on the panel's own
# min/max (unlike scales::breaks_pretty(n = 3), which anchors near 0 and
# will happily place a break at or past a data extreme when the range
# doesn't start near zero -- the X11 dual-axis figures' actual complaint).
# Searches the largest "nice" step (a single leading digit x 10^k, e.g. 300,
# 2000 -- deliberately broader than nice_step_for_target()'s 1/2/2.5/5
# family, since that family alone can't always reach 3 ticks that both fit
# inside the range AND use as much of it as this does) for which 3 evenly
# spaced multiples of it fit strictly inside (lo, hi); falls back to one
# order of magnitude down (repeatedly) if nothing fits, since a data range
# with very little "room" (e.g. lo and hi close to a shared magnitude's
# multiples) can otherwise have no valid step at the initial magnitude.
# Returns a breaks-function, like scales::breaks_pretty(), so it drops
# straight into scale_y_continuous(breaks = ...)/sec_axis(breaks = ...).
three_ticks_away_from_limits <- function(){
  function(limits){
    lo <- limits[1]; hi <- limits[2]; span <- hi - lo
    if (!is.finite(span) || span <= 0) return(scales::breaks_pretty(n = 3)(limits))

    magnitude <- 10 ^ floor(log10(span / 3))
    for (i in 0:4) {
      mag <- magnitude / 10^i
      for (step in mag * 9:1) {
        lowest <- ceiling(lo / step + 1e-9) * step
        ticks <- lowest + step * 0:2
        if (all(ticks > lo) && all(ticks < hi)) return(ticks)
      }
    }
    scales::breaks_pretty(n = 3)(limits) # fallback, shouldn't normally trigger
  }
}

# Natural Earth 10m "Land" (naturalearthdata.com, public domain), vendored
# into data/EUROPE_shapefile/ 2026-08-11: a single global land-polygon layer
# at the same 10m-class resolution as the GADM France file this replaced,
# covering all of Europe and the UK (and everywhere else) in one file. Was
# previously France-only (GADM ADM0), overplotted on the crude `maps`-
# package "world" coastline for every other country -- visibly mismatched
# resolution at the UK/Belgium/Spain edges of every map. No per-country
# subsetting needed: "Land" is one continuous landmass polygon, matching
# how it's used below (a single filled black layer, not coloured by
# country), so create_the_basic_map() no longer needs the low-res
# map_data("world") layer at all when high_res_coast = TRUE.
.high_res_coastline_cache <- NULL
high_res_coastline <- function(){
  if(is.null(.high_res_coastline_cache)){
    library(sf)
    .high_res_coastline_cache <<- sf::st_read("data/EUROPE_shapefile/ne_10m_land.shp", quiet = TRUE) |>
      sf::st_coordinates() |>
      as.data.frame() |>
      dplyr::transmute(long = X, lat = Y, group = interaction(L1, L2, L3, drop = TRUE))
  }
  .high_res_coastline_cache
}

create_the_basic_map <- function(map_df, var_name,
                                 in_situ_fixed_station = NULL,
                                 cruise_stations = NULL,
                                 glider_stations = NULL,
                                 legend_limits = NULL,
                                 log_scale = TRUE,
                                 high_res_coast = FALSE) {
  
  if (str_detect(var_name, 'chl|CHL')) {
    title = "Chl-a"
    unit = "mg m⁻³"
    if (legend_limits |> is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  if (str_detect(var_name, 'tsm|SPM|TSM|plume')) {
    title = "SPM"
    unit = "g m⁻³"
    # legend_limits <- map_df$analysed_spim[which(map_df$plume)] |> quantile(probs = c(0.1, 0.9), na.rm = TRUE)
    if (legend_limits |> is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  FRANCE_shapefile <- map_data('world')[map_data('world')$region == "France",]
  
  the_base_map <- ggplot() + 
    geom_raster(data = map_df, aes(x = lon, y = lat, fill = analysed_spim), interpolate = FALSE) + 
    scale_fill_viridis_c(na.value = "transparent", option = "viridis", trans = if(log_scale) "log10" else "identity",
                         limits = c(legend_limits[1], legend_limits[2]), oob = scales::squish,
                         n.breaks = 5, name = paste(title, " (", unit, ")", sep = "")) +
    guides(fill = guide_colourbar(title.position = "right"))
  
  if (var_name == 'plume') {
    the_base_map <- the_base_map + geom_raster(data = map_df[which(map_df$plume),], aes(x = lon, y = lat), fill = "red", interpolate = FALSE) 
  }
  
  the_map <- the_base_map +

    (if(high_res_coast){
      # High-resolution coastline (Natural Earth 10m land, all of Europe +
      # UK and beyond -- see high_res_coastline() above) as the only layer.
      # Fixed 2026-08-11: previously this was France-only (GADM), layered
      # OVER the crude low-res map_data("world") coastline underneath, so
      # every neighbouring country's coast (UK, Belgium, Spain) stayed at
      # visibly cruder resolution than France's own -- both on the main
      # national panel and every zoomed inset. The new file covers the
      # whole extent any of these maps needs, so the low-res layer is no
      # longer drawn at all when high_res_coast = TRUE.
      list(
        geom_polygon(data = high_res_coastline(), aes(x = long, y = lat, group = group), color = 'grey60', fill = 'black')
      )
    } else {
      list(
        # First layer: worldwide map
        geom_polygon(data = map_data("world"), aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black'),
        # Second layer: Country map
        geom_polygon(data = FRANCE_shapefile, aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black')
      )
    }) +
    coord_cartesian(xlim = range(map_df$lon), ylim = range(map_df$lat), expand = FALSE) +
    
    scale_x_continuous(name = "", labels = function(x) paste(x, "°E", sep = "")) +
    scale_y_continuous(name = "", labels = function(x) paste(x, "°N", sep = "")) +
    ggplot_theme() + 
    
    theme(plot.title = element_text(size = 45),
          legend.position = "right",
          legend.title = element_text(angle = -90, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key.height = unit(6, "lines"),
          legend.key.width = unit(3, "lines")) 
  
  return(the_map)
  
}

plot_x11_river_and_plume <- function(X11_data, type_of_signal, show_axis_titles = TRUE) {

  unique_years <- X11_data$dates |> year() |> unique()

  X11_data_for_plot <- X11_data |>
    rename(river_flow = !!sym(paste(type_of_signal, "signal_river_flow", sep = "_")),
           plume_area = !!sym(paste(type_of_signal, "signal_plume_area", sep = "_"))) |>
    select(dates, river_flow, plume_area)

  if (type_of_signal %in% c("Seasonal", "Residual")) {
    X11_data_for_plot <- X11_data_for_plot |>
      mutate(river_flow = river_flow + mean(X11_data$Raw_signal_river_flow, na.rm = T),
             plume_area = plume_area + mean(X11_data$Raw_signal_plume_area, na.rm = T))
  }

  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = X11_data_for_plot$river_flow,
                                                 var_ref = X11_data_for_plot$plume_area)

  X11_data_for_plot <- X11_data_for_plot |> mutate(river_flow_scaled = river_flow * scaling_factor$diff + scaling_factor$adjust)

  r_value <- cor(X11_data_for_plot$plume_area, X11_data_for_plot$river_flow, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  the_plot <- ggplot() +

    geom_path(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown", linewidth = 1) +

    geom_path(data = X11_data_for_plot, aes(x = dates, y = river_flow_scaled), color = "blue", linewidth = 1) +

    annotate("text", x = min(X11_data_for_plot$dates), y = Inf, label = r_label,
            hjust = 0, vjust = 1.5, size = 6, colour = "black") +

    scale_x_date(name = "",
                 breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                 labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist()) +

    scale_y_continuous(name = if(show_axis_titles) "Plume area (km²)" else NULL,
                       breaks = three_ticks_away_from_limits(),
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = if(show_axis_titles) "River flow (m³ s⁻¹)" else NULL,
                                           breaks = three_ticks_away_from_limits())) +

    labs(title = paste(type_of_signal, "signal")) +
    ggplot_theme() +
    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(size=30, colour = "black"),
          text = element_text(size=25, colour = "black"),
          axis.text = element_text(size=20, colour = "black"),
          axis.title = element_text(size=30, colour = "black")) +
    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y.left = element_text(color = "brown"), 
          axis.ticks.y.left = element_line(color = "brown"),
          axis.line.y.left = element_line(color = "brown"),
          axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
          
          axis.text.y.right = element_text(color = "blue"), 
          axis.ticks.y.right = element_line(color = "blue"),
          axis.line.y.right = element_line(color = "blue"),
          axis.title.y.right = element_text(color = "blue", margin = unit(c(0, 0, 0, 7.5), "mm")),
          
          panel.border = element_rect(linetype = "solid", fill = NA))
  
  return(the_plot)
  
}


