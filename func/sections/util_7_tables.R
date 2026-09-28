# Tables ---------------------------------------------------------------

# Create pretty tables
# file_path = "output/MATCH_UP_DATA/FRANCE/STATISTICS/"; sat_name = "SEXTANT"
validation_tables <- function(file_path, sat_name) {
  
  # Load all output stats
  # TODO: Clean up zone names and order by latitude
  stat_files <- map_dfr(dir(file_path, pattern = paste0(sat_name,"_stats_all"), full.names = TRUE), read_csv) |> 
    filter(zone != "GLOBAL", source == "ALL", site == "ALL", season == "ALL") |> 
    dplyr::select(zone, Bias, Error, n, variable) |> 
    mutate_at(vars(-zone, -variable), ~ round(., 1)) |> 
    mutate(var_label = paste("When compared with in situ", ifelse(variable == 'TUR', 'TURB', variable))) |>
    dplyr::select(-variable)
  # dplyr::rename(`Bias (%)` = Bias, `Error (%)` = Error)
  
  desired_colnames <- names(stat_files) |> 
    str_replace_all("Bias", "Bias (%)") |> 
    str_replace_all("Error", "Error (%)") |> 
    str_remove_all("var_label")
  names(desired_colnames) <- names(stat_files)
  
  # Create the SPM/TUR table
  # TODO: Wrap this into a function to be called per variable
  table_SPM_TUR <- stat_files |>
    filter(!grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite SPM', subtitle = 'Compared with SPM and Turbidity in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    # Pretty tabs
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3), # A bit more padding for summaries
                row_group.padding = px(4)) |> # And even more for our groups
    # More styling
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Create the CHLA table
  table_CHLA <- stat_files |>
    filter(grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite CHL', subtitle = 'Compared with CHLA in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3),
                row_group.padding = px(4)) |> 
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Save and exit
  ## SPM/TUR
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.html"), inline_css = TRUE)
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.png"), expand = 10)
  ## CHLA
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.html"), inline_css = TRUE)
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.png"), expand = 10)
}

