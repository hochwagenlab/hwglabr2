#' Return elapsed time
#'
#' Given a time stamp generated with \code{proc.time}, returns elapsed time in
#' sec., min. or hrs, as appropriate.
#' @param t0 Output of \code{proc.time}, either all elements or just element
#' "elapsed". No default.
#' @return \code{String}.
#' @examples
#' \dontrun{
#' start <- proc.time()
#' elapsed_time(t0=start)
#' 
#' start <- proc.time()[3]
#' elapsed_time(t0=start)
#' }
#' @export

elapsed_time <- function(t0){
  
  # IO checks
  if (!(is(t0, "proc_time"))) {
    if (!(is(t0, "numeric") | length(t0) != 1)) {
      stop("'t0' must be the output of proc.time (in full or a single element)")
    }
  }
  
  # Are all elements returned by proc.time included?
  if (is(t0, "proc_time")) t0 <- t0[3]
  # Calculate elapsed time
  elapsed_time <- proc.time()[3] - t0
  
  # Convert to appropriate unit
  if (elapsed_time < 60) {
    elapsed_time <- paste0(round(elapsed_time, 1), " sec.")
  } else if (elapsed_time >= 60 & elapsed_time < 3600) {
    elapsed_time <- paste0(round(elapsed_time / 60, 1), " min.")
  } else {
    elapsed_time <- paste0(round(elapsed_time / 60 / 60, 1), " hrs.")
  }
  return(elapsed_time)
}
