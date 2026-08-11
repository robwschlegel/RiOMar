# One-off: generate the LaTeX longtable BODY for the new Supplementary
# "compass-octant occupancy trend" table (\label{tab:s_octant_trends} in
# manuscript.tex), from output/STATS/octant_trend_compact_summary.csv
# (func/compute_direction_octant_trend.R). Mirrors
# generate_table_s_monthly_trends.R's fmt_cell()/precision/significance-star
# convention and Table 7's (\label{tab:monthly_trends}) "sig. months (dir.)"
# row, here with an added row for the annual-level trend, one block of 2
# rows per (driver, octant). Unlike Table 7 (a page-fitting table wrapped in
# \resizebox), this is a multi-page longtable, which can't be resized the
# same way -- the monthly min/max range that Table 7 also prints is
# deliberately omitted here (full detail stays in the CSV only) because
# printing it caused a genuine page-width overflow (four "min to max"
# strings side by side), not merely a cosmetic one.
#
# Slopes in octant_trend_compact_summary.csv are proportion-of-days-per-day
# (from lm(proportion ~ date)); converted here to percentage points per year
# (pp yr^-1) for readability, matching the "annualise the per-day slope"
# convention already used for every other trend table in this pipeline.
#
# Run from repo root: Rscript func/generate_table_s_octant_trends.R
suppressMessages(library(tidyverse))

YR <- 365.25
precision <- 2  # pp yr^-1

drivers <- tibble::tribble(
  ~driver,    ~label,
  "wind",     "Wind",
  "wave",     "Wave",
  "current",  "Current"
)
octants <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
zones <- c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION")

fmt_cell <- function(slope, p, precision){
  if(is.na(slope)) return("--")
  sign <- if(slope >= 0) "$+$" else "$-$"
  val <- sprintf(paste0("%.", precision, "f"), abs(slope))
  star <- if(!is.na(p) && p < 0.05) "$^{*}$" else ""
  paste0(sign, val, star)
}

detail <- read_csv("output/STATS/octant_trend_compact_summary.csv", show_col_types = FALSE)

lines <- c()
for(i in seq_len(nrow(drivers))){
  drv <- drivers$driver[i]; label <- drivers$label[i]

  for(oct in octants){
    row_annual <- purrr::map_chr(zones, function(z){
      r <- dplyr::filter(detail, driver == drv, zone == z, octant == oct)
      if(nrow(r) == 0 || is.na(r$annual_slope)) return("--")
      fmt_cell(r$annual_slope * YR * 100, r$annual_slope_p, precision)
    })
    row_sig <- purrr::map_chr(zones, function(z){
      r <- dplyr::filter(detail, driver == drv, zone == z, octant == oct)
      if(nrow(r) == 0 || is.na(r$n_significant_months)) return("--")
      # "mixed" abbreviated to "mix" here (print layer only) purely for
      # longtable column width -- octant_trend_compact_summary.csv itself
      # keeps the full "mixed" string.
      dir_label <- if(r$monthly_majority_direction == "mixed") "mix" else r$monthly_majority_direction
      paste0(r$n_significant_months, "/", r$n_months_fit, " (", dir_label, ")")
    })

    lead <- paste0(label, " (", oct, ")")
    lines <- c(lines,
               paste0(lead, " & annual & ", paste(row_annual, collapse = " & "), " \\\\"),
               paste0(" & sig./12 (dir.) & ", paste(row_sig, collapse = " & "), " \\\\"))
  }
  if(i < nrow(drivers)) lines <- c(lines, "\\addlinespace")
}

out_path <- "output/STATS/table_s_octant_trends_body.tex"
writeLines(lines, out_path)
cat("Wrote", length(lines), "lines to", out_path, "\n")
cat("Paste between \\endhead and \\end{longtable} in manuscript.tex (\\label{tab:s_octant_trends}).\n")
