# func/config.R
# Project-wide settings, read from metadata/riomar_config.yml. The Python twin
# is func/config.py -- keep the data-root resolution logic identical in both.


# Setup ------------------------------------------------------------------

riomar_config_path <- "metadata/riomar_config.yml"


# Functions ---------------------------------------------------------------

riomar_config <- function(){
  yaml::read_yaml(riomar_config_path)
}

# Absolute path of the large-dataset root (SEXTANT, WIND, WAVE, GLORYS...).
# Resolution order: RIOMAR_DATA_ROOT env var > `data_root` in the YAML > the
# OS's pCloud folder (~/pCloud Drive/data on macOS, ~/pCloudDrive/data
# elsewhere), falling back to the other spelling if the OS default is missing.
riomar_data_root <- function(){
  explicit <- Sys.getenv("RIOMAR_DATA_ROOT", unset = "")
  if (explicit == "") explicit <- riomar_config()$data_root
  if (!is.null(explicit) && explicit != "") {
    candidates <- explicit
  } else {
    pcloud_folders <- c("pCloud Drive", "pCloudDrive")
    if (Sys.info()[["sysname"]] != "Darwin") pcloud_folders <- rev(pcloud_folders)
    candidates <- file.path("~", pcloud_folders, "data")
  }
  for (candidate in candidates) {
    path <- normalizePath(path.expand(candidate), mustWork = FALSE)
    if (dir.exists(path)) return(path)
  }
  stop("RiOMar data root not found; tried: ", paste(candidates, collapse = ", "),
       ". Set RIOMAR_DATA_ROOT or data_root in metadata/riomar_config.yml.")
}

# riomar_data_path("WIND", zone) -> <data_root>/WIND/<zone>
riomar_data_path <- function(...){
  file.path(riomar_data_root(), ...)
}

riomar_zones <- function(){
  unlist(riomar_config()$zones)
}
