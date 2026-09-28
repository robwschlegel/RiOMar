# Statistics --------------------------------------------------------------

# Check for leap year
is_leap_year <- function(year) {
  (year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)
}

# Function to adjust DOY for non-leap years
adjust_doy <- function(year, doy) {
  if (!is_leap_year(year) && doy > 59) {
    doy + 1
  } else {
    doy
  }
}

# Lagged correlations
# NB: The x vector is lagged, meaning it is effectively pushed forward in time
lagged_correlation <- function(x, y, max_lag) {
  df_cor <- tibble(
    lag = 0:max_lag,
    cor = map_dbl(0:max_lag, ~ cor(y, lag(x, .), use = "pairwise.complete.obs"))
  )
  return(df_cor)
}

# Calculate STL
stl_single <- function(x_col, out_col, start_date, ts_freq = 365.25){
  
  # Create ts object and calculate stl
  ts_x <- ts(zoo::na.approx(x_col), frequency = ts_freq, start = c(year(start_date), quarter(start_date)))
  stl_x <- stl(ts_x, s.window = "periodic", t.window = 11, inner = 5, outer = 0)
  #, robust = TRUE) # rather allow outliers to influence results
  
  # Add NA to end if necessary
  if(length(stl_x$time.series[,2]) != length(x_col)){

    # Create a matrix of NA values
    na_matrix <- matrix(NA, nrow = length(x_col)-length(stl_x$time.series[,2]), ncol = ncol(stl_x$time.series))
    
    # Append the NA matrix to the original matrix
    stl_x$time.series <- rbind(stl_x$time.series, na_matrix)
  }
  
  # Decide which column to take
  if(out_col == "seas"){
    return(as.vector(stl_x$time.series[,1]))
  } else if(out_col == "inter"){
    return(as.vector(stl_x$time.series[,2]))
  } else if(out_col == "remain"){
    return(as.vector(stl_x$time.series[,3]))
  } else {
    stop("'out_col' not recognised")
  }
}

# The full suite of stats to calculate
compute_stats <- function(x_vec, y_vec){
  
  if(!is.numeric(x_vec)) stop("x_vec is not numeric")
  if(!is.numeric(y_vec)) stop("y_vec is not numeric")
  
  if(length(x_vec) < 3){
    return(data.frame(row.names = NULL,
                      n = length(x_vec),
                      Slope = NA,
                      Slope_log = NA,
                      RMSE = NA,
                      MSA = NA,
                      MAPE = NA,
                      Bias = NA,
                      Error = NA))
  }
  
  # Calculate RMSE (Root Mean Square Error)
  rmse <- sqrt(mean((y_vec - x_vec)^2, na.rm = TRUE))
  
  # Calculate MAPE (Mean Absolute Percentage Error)
  mape <- mean(abs((y_vec - x_vec) / x_vec), na.rm = TRUE) * 100
  
  # Calculate MSA (Mean Squared Adjustment)
  msa <- mean(abs(y_vec - x_vec), na.rm = TRUE)
  
  # Calculate linear slope
  lin_fit <- lm(y_vec ~ x_vec)
  slope <- coef(lin_fit)[2]
  
  # Calculate log-log linear slope
  log_lin_fit <- lm(log10(y_vec) ~ log10(x_vec))
  log_slope <- coef(log_lin_fit)[2]
  
  # Calculate Bias
  log_ratio <- log10(y_vec / x_vec)
  log_ratio_median <- median(log_ratio, na.rm = TRUE)
  bias_perc <- 100 * (sign(log_ratio_median) * (10^abs(log_ratio_median) - 1))
  
  # Calculate error
  log_ratio_median_abs <- median(abs(log_ratio), na.rm = TRUE)
  error_perc <- 100 * (10^log_ratio_median_abs - 1)
  
  # Combine int data.frame and exit
  return(data.frame(row.names = NULL,
                    n = length(x_vec),
                    Slope = round(slope, 2),
                    Slope_log = round(log_slope, 2),
                    RMSE = round(rmse, 6),
                    MSA = round(msa, 6),
                    MAPE = round(mape, 2),
                    Bias = round(bias_perc, 2),
                    Error = round(error_perc, 2)))
}


