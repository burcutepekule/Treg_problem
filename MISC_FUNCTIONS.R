get_last_below_threshold_run = function(x, threshold) {
  below <- x < threshold
  rle_result <- rle(below)
  lengths <- rle_result$lengths
  values  <- rle_result$values
  ends    <- cumsum(lengths)
  
  # Find runs where value is TRUE (i.e. below threshold)
  last_run_index <- max(which(values == TRUE))
  
  # Get start and end indices of the last run
  end_index <- ends[last_run_index]
  start_index <- end_index - lengths[last_run_index] + 1
  
  # return(list(indices = start_index:end_index, values = x[start_index:end_index]))
  return(start_index)
}

get_steady_time =  function(x, smoothing_window, roll_window, threshold) {
  smoothed = zoo::rollmean(x, k = smoothing_window, fill = NA, align = "center")
  rolling_sd = RcppRoll::roll_sd(smoothed, n = roll_window, fill = NA, align = "right")
  index_steady = get_last_below_threshold_run(rolling_sd, threshold = threshold)
  return(index_steady)
}