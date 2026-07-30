# func/surface.R
# Suspicious plume-day flagging
#
# Scans every panache Results.csv (both dynamic- and static-threshold runs, all four zones) 
# for single/few-day spikes in plume area that are almost certainly detection artefacts 
# rather than real events.
#
# Motivation: the combined plume-area time series (Figure_10-style plot) shows
# a handful of days per zone where area jumps to several times the seasonal
# maximum for a single day and then collapses back to normal. Visual
# inspection of the daily plume-shape grids (one panel per month x year)
# confirms these are panache's flood-fill briefly spilling out past the river
# mouth into open-ocean pixels -- most plausible when river-mouth SPM is only
# barely above ambient/open-ocean SPM that day, so the (dynamic or static)
# threshold no longer discriminates plume from background. This is
# distinct from a genuine flood/storm event, which builds and decays over
# several days and therefore never produces the extreme single-day
# spike-to-neighbourhood ratio used below to flag it.

source("func/multi.R")


# Loading -------------------------------------------------------------------

# Read one zone/threshold Results.csv, gap-filled to a continuous daily
# sequence. Keeps every core column (date through confidence_index_in_perc,
# same slice used by util.R::load_plume_ts()) but, unlike load_plume_ts(),
# does NOT null out large values so nothing gets silently discarded before 
# it can be flagged. The per-river SPM_threshold_* columns are dropped since their
# names differ by zone (one column per river mouth in that zone).
# load_plume_results("GULF_OF_LION", "dynamic")
load_plume_results <- function(zone, threshold = c("dynamic", "static")){
  threshold <- match.arg(threshold)
  file_name <- paste0("output/panache/", threshold, "/", zone, "/Results.csv")
  suppressMessages({
    df <- read_csv(file_name) |>
      dplyr::mutate(date = as.Date(date)) |>
      dplyr::select(date:confidence_index_in_perc) |>
      tidyr::complete(date = seq(min(date), max(date), by = "day"))
  })
  df |> dplyr::mutate(zone = zone, threshold = threshold, .before = "date")
}


# Flagging --------------------------------------------------------------------

# Flag a day as suspicious if its plume area is both (a) in the extreme tail
# of that zone/threshold's own area distribution (area > the p_floor
# quantile) and (b) many times larger than the surrounding +/-window_days
# neighbourhood's median (a leave-one-out median, so the flagged day itself
# never inflates its own baseline). Both conditions matter: (a) alone flags
# genuine multi-day flood peaks (which sit in the tail too, just not
# isolated); (b) alone flags trivial blips in otherwise near-zero, low-flow
# stretches, where any small bump is "many times" a tiny local median.
#
# ratio_min = 15 was chosen by hand against the two events found by visual
# inspection of the daily plume-shape grids -- BAY_OF_SEINE 2014-02-02
# (9563 km^2, ratio 61 against its neighbourhood) and GULF_OF_LION
# 2005-05-12 (22359 km^2, ratio 203) -- versus the largest genuine flood
# events in the 1998-2025 record (e.g. GULF_OF_LION 2018-03-01/02,
# BAY_OF_BISCAY 2018-01-15/17), whose ratios top out below 16. 15 sits in
# the gap between those two populations across all four zones and both
# threshold runs.
# flag_suspicious_plume_days("GULF_OF_LION", "dynamic")
flag_suspicious_plume_days <- function(zone, threshold = c("dynamic", "static"),
                                        window_days = 7, ratio_min = 15, p_floor = 0.95){
  threshold <- match.arg(threshold)
  df <- load_plume_results(zone, threshold)
  area <- df$area_of_the_plume_mask_in_km2
  n <- length(area)

  local_median <- vapply(seq_len(n), function(i){
    idx <- setdiff(max(1, i - window_days):min(n, i + window_days), i)
    stats::median(area[idx], na.rm = TRUE)
  }, numeric(1))

  floor_val <- stats::quantile(area, p_floor, na.rm = TRUE)
  ratio <- area / pmax(local_median, 1)

  df |>
    dplyr::mutate(local_median = local_median, ratio = ratio,
                  suspicious = !is.na(area) & area > floor_val & ratio > ratio_min)
}

# Run flag_suspicious_plume_days() over every zone and both threshold runs.
# flag_all_suspicious_plume_days()
flag_all_suspicious_plume_days <- function(window_days = 7, ratio_min = 15, p_floor = 0.95){
  purrr::map_dfr(zones, function(zone){
    purrr::map_dfr(c("dynamic", "static"), function(threshold){
      flag_suspicious_plume_days(zone, threshold, window_days, ratio_min, p_floor)
    })
  })
}

# Plotting --------------------------------------------------------------------

# Daily plume-area time series for all four zones, with suspicious days
# overplotted in red. One plot per threshold run (dynamic areas and static
# areas for the same zone are not on a comparable scale, so they don't share
# a panel). Mirrors the daily time-series layout in
# multi.R::multi_plot()'s plot_daily() -- facet_wrap(~zone, ncol = 1,
# scales = "free_y"), one panel per zone, year-labelled x-axis.
# plot_suspicious_plume_days(flag_all_suspicious_plume_days(), "dynamic")
plot_suspicious_plume_days <- function(df_flagged, threshold = c("dynamic", "static")){
  threshold <- match.arg(threshold)
  df <- dplyr::filter(df_flagged, threshold == !!threshold)
  unique_years <- df$date |> lubridate::year() |> unique()

  pl <- ggplot(df, aes(x = date, y = area_of_the_plume_mask_in_km2)) +
    geom_path(color = "grey40") +
    geom_point(color = "grey40", size = 0.8) +
    geom_point(data = dplyr::filter(df, suspicious), color = "red", size = 2) +
    facet_wrap(~zone, ncol = 1, scales = "free_y") +
    scale_x_date(name = "", expand = c(0.01, 0.01),
                 breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                 labels = unique_years |> stringr::str_extract_all('[0-9][0-9]$') |> unlist()) +
    scale_y_continuous(name = "Plume area (km²)") +
    labs(title = paste0(stringr::str_to_title(threshold), " threshold")) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  ggsave(filename = paste0("figures/qc/suspicious_plume_days_", threshold, ".png"),
        plot = pl, width = 24, height = 24, dpi = 300)
  invisible(pl)
}

# Run plot_suspicious_plume_days() for both threshold runs off one shared
# flag_all_suspicious_plume_days() call.
# plot_all_suspicious_plume_days()
plot_all_suspicious_plume_days <- function(window_days = 7, ratio_min = 15, p_floor = 0.95){
  df_flagged <- flag_all_suspicious_plume_days(window_days, ratio_min, p_floor)
  plot_suspicious_plume_days(df_flagged, "dynamic")
  plot_suspicious_plume_days(df_flagged, "static")
  invisible(df_flagged)
}

# Run it
flag_all_suspicious_plume_days()
