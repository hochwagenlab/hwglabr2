#' Return elapsed time
#'
#' Given a start and end time in seconds, returns elapsed time in sec., min. or
#' hrs, as appropriate.
#' @param t0 Numeric start time representing seconds. No default.
#' @param t1 Numeric end time representing continuous second count from
#' \code{t0}. No default.
#' @return \code{String} containing elapsed time.
#' @examples
#' \dontrun{
#' start <- proc.time()[3]
#' end <- proc.time()[3]
#' elapsed_time(t0=start, t1=end)
#' 
#' start <- proc.time()
#' elapsed_time(start[3], proc.time()[3])
#' }
#' @export

elapsed_time <- function(t0, t1){
  
  # IO checks
  if (!(is(t0, "numeric") && is(t1, "numeric"))) {
    stop('"t0" and "t1" must be numeric')
  }
  
  # Calculate elapsed time
  elapsed_time <- t1 - t0
  
  # Convert to appropriate unit
  if (elapsed_time < 60) {
    elapsed_time <- paste0(round(elapsed_time, 1), " sec.")
  } else if (elapsed_time >= 60 & elapsed_time < 3600) {
    elapsed_time <- paste0(round(elapsed_time / 60, 1), " min.")
  } else {
    elapsed_time <- paste0(round(elapsed_time / 60 / 60, 2), " hrs.")
  }
  return(elapsed_time)
}
