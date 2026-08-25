# One-off: generate the LaTeX longtable BODY for the Supplementary "full
# month-by-month linear trend detail" table (\label{tab:s_monthly_trends} in
# manuscript.tex), from output/STATS/monthly_trend_summary.csv
# (func/compute_seasonal_trend.R). Dynamic threshold only, AR(1)/HAC-weighted
# ("ar") trend -- matching the compact main-text Table 6
# (\label{tab:monthly_trends}) convention this table is the full detail
# behind.
#
# This project hand-transcribes every table into manuscript.tex rather than
# \input{}-ing a generated fragment (see manuscript/TODO.md's "Reference:
# where each table's numbers come from" note) -- this script follows that
# convention too: it writes a reviewable .tex fragment, not a live \input.
# Paste the fragment's contents into manuscript.tex between the longtable's
# \endhead and \end{longtable}, replacing the existing hand-transcribed rows.
#
# Run from repo root: Rscript func/generate_table_s_monthly_trends.R
suppressMessages(library(tidyverse))

YR <- 365.25

# variable key -> (row label as it appears in manuscript.tex, decimal
# precision, unit conversion factor applied to the annualised slope).
# SPM_mass is stored in kg (see func/compute_mass_spm_trend.R's comment);
# manuscript tables report it in tonnes.
variables <- tibble::tribble(
  ~variable,        ~label,                                  ~precision, ~unit_convert,
  "plume_area",     "Plume area (km$^2$ yr$^{-1}$)",          1,          1,
  "SPM_mass",       "SPM mass (t yr$^{-1}$)",                 1,          1 / 1000,
  "compactness",    "Compactness (yr$^{-1}$)",                4,          1,
  "alongcoast_km",  "Along-coast (km yr$^{-1}$)",             3,          1,
  "flow",           "Discharge (m$^3$ s$^{-1}$ yr$^{-1}$)",   1,          1,
  "wind",           "Wind (m s$^{-1}$ yr$^{-1}$)",            3,          1,
  "tide",           "Tide (m yr$^{-1}$)",                     4,          1,
  "wave",           "Wave (m yr$^{-1}$)",                     4,          1,
  "current",        "Current (m s$^{-1}$ yr$^{-1}$)",         4,          1
)

# Column order matches the manuscript table's zone columns left to right.
zones <- c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION")

fmt_cell <- function(slope, p, precision){
  if(is.na(slope)) return("--")
  sign <- if(slope >= 0) "$+$" else "$-$"
  val <- sprintf(paste0("%.", precision, "f"), abs(slope))
  star <- if(!is.na(p) && p < 0.05) "$^{*}$" else ""
  paste0(sign, val, star)
}

detail <- read_csv("output/STATS/monthly_trend_summary.csv", show_col_types = FALSE) |>
  dplyr::filter(threshold == "dynamic", weight_choice == "ar")

lines <- c()
for(i in seq_len(nrow(variables))){
  v <- variables$variable[i]; label <- variables$label[i]
  precision <- variables$precision[i]; conv <- variables$unit_convert[i]

  for(m in 1:12){
    row_vals <- purrr::map_chr(zones, function(z){
      r <- dplyr::filter(detail, variable == v, zone == z, month == m)
      if(nrow(r) == 0 || is.na(r$slope)) return("--")
      fmt_cell(r$slope * YR * conv, r$slope_p, precision)
    })
    lead <- if(m == 1) label else ""
    lines <- c(lines, paste0(lead, " & ", month.abb[m], " & ", paste(row_vals, collapse = " & "), " \\\\"))
  }
  if(i < nrow(variables)) lines <- c(lines, "\\addlinespace")
}

out_path <- "output/STATS/table_s_monthly_trends_body.tex"
writeLines(lines, out_path)
cat("Wrote", length(lines), "lines to", out_path, "\n")
cat("Paste between \\endhead and \\end{longtable} in manuscript.tex (\\label{tab:s_monthly_trends}), replacing the existing hand-transcribed rows.\n")
