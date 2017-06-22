#' Mean, SD, SEM, and confidence interval
#'
#' Calculates the mean value and confidence interval, as well as standard
#' deviation and standard error of the mean, for every column in an input table.
#' Confidence intervals are calculated using the "Empirical bootstrap" technique
#' (based on random sampling with replacement). The input will typically be a
#' table of signal values for different obervations (across rows) over a set of
#' positions (across columns).
#' 
#' @param signal_data Input signal data matrix (or an object that can be coerced
#' to class \code{matrix}). No default.
#' @param ci Numeric value between 0 and 1 representing the confidence interval
#' percentage. Defaults to \code{0.95}.
#' @param rep_bootstrap Integer specifying the number of times to repeat the
#' bootstrapping procedure. Defaults to \code{1000}.
#' @param na_rm Logical indicating whether to remove \code{NA} values in each
#' column before calculating the mean and bootstrapping. Defaults to \code{TRUE}.
#' @return A matrix with as many rows as there are columns in the input and
#' three columns:
#' \enumerate{
#'   \item \code{Mean} Input column mean
#'   \item \code{Lower} Lower limit of the bootstrapping confidence interval
#'   \item \code{Upper} Upper limit of the bootstrapping confidence interval
#'   \item \code{SD} Standard deviation
#'   \item \code{SEM} Standard error of the mean
#' }
#' @examples
#' \dontrun{
#' signal_mean_and_ci(WT_signal_at_orf)
#' 
#' signal_mean_and_ci(WT_signal_at_summits, ci=0.90)
#' 
#' signal_mean_and_ci(WT_signal_at_dsb, bootstrap_reps=5000)
#' }
#' @export

signal_mean_and_ci <- function(signal_data, ci=0.95, rep_bootstrap=1000,
                               na_rm=TRUE) {
  t0  <- proc.time()[3]
  
  # IO checks
  signal_data <- as.matrix(signal_data)
  
  if (!is(signal_data[1, ], "numeric")) {
    stop('"signal_data" must be either a numeric matrix or',
         'an object that can be coerced to one', call. = FALSE)
  }
  
  message('Bootstrapping ', rep_bootstrap, ' times...')
  message('(confidence interval = ', ci, ')')
  
  # Apply function mean_and_ci (defined here below) to every column in the input
  result <- apply(signal_data, 2, bootstrap_mean_ci, method='empirical',
                  conf_int=ci, B=rep_bootstrap,
                  na.rm=na_rm)
  
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
  return(t(result))
}

bootstrap_mean_ci <- function(x, method='empirical', conf_int=0.95,
                              B=1000, na.rm=TRUE) {
  stopifnot(method %in% c('percentile', 'pivot', 'empirical'))
  
  if(na.rm) x <- x[!is.na(x)]
  
  n <- length(x)
  theta <- mean(x)
  std <- sd(x)
  sem <- std / sqrt(n)
  
  if (n < 2) return(c(Mean=xbar, Lower=NA, Upper=NA))
  
  # Draw random samples with replacement the same size as the input
  # as many times as specified (the x-star samples) and compute their means theta-hat
  #z <- replicate(B, mean(sample(x, n, replace=TRUE), na.rm = T))
  # More efficient implementation (from the "Hmisc" package)
  theta_hat <- unlist(lapply(seq_len(B),
                            function(i, x, N) sum(x[sample.int(n=N, size=N,
                                                               replace=TRUE)]),
                            x=x, N=n)) / n
  
  percentiles <- quantile(theta_hat, c((1 - conf_int) / 2, (1 + conf_int) / 2))
  
  if (method == 'percentile') {
    limits <- percentiles
  } else if (method == 'pivot') {
    limits <- rev(2 * theta - percentiles)
  } else if (method == 'empirical') {
    deltas <- theta_hat - theta
    percentiles <- quantile(deltas, c((1 - conf_int) / 2, (1 + conf_int) / 2))
    limits <- theta - rev(percentiles)
  }
  
  names(limits) <- NULL
  return(c(Mean=theta, Lower=limits[1], Upper=limits[2], SD=std, SEM=sem))
}
