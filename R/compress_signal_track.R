#' Average signal variable over specified position windows
#'
#' Partitions the position variable (first column of input data frame) in
#' windows of the size specified in \code{window_size} and computes the average
#' signal (second column) in each partition. Most common use case is to compress
#' the data for plotting.
#' @param df Input signal track data as a \code{data.frame} with at least two
#' columns:
#' \enumerate{
#'   \item Numeric column with positions to partition
#'   \item Numeric column with values to average in each partition
#' }
#' No default.
#' @param window_size Integer indicating the length (in bp) of the window (or
#' partition) to average signal in. No default.
#' @examples
#' \dontrun{
#' compress_signal_track(df=signal_data, window_size=10)
#' 
#' compress_signal_track(chrI_signal, 100)
#' }
#' @export

compress_signal_track <- function(df, window_size) {
  t0 <- proc.time()[3]
  
  # IO checks
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' required for this function to work.\n", 
         "Please install it: install.packages('dplyr')", call. = FALSE)
  }
  
  colnames(df)[1:2] <- c('position', 'signal')
  #The two first columns are used; must be numeric
  if (!is(df$position, "numeric") | !is(df$signal, "numeric")) {
    stop('The two first columns of the input object must be numeric.',
         call. = FALSE)
  }

  message('Build tiles of the position variable...')
  df$tile <- make_virtual_tiles(df$position, window_size)
  
  message('Compute mean per tile...')
  df <- dplyr::summarise(dplyr::group_by(df, tile),
                         position=tile_midpoint(position),
                         window_mean=mean(signal, na.rm = T))
  
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
  return(as.data.frame(df[, c('position', 'window_mean')]))
}

tile_midpoint <- function(x) mean(c(min(x), max(x)))

make_virtual_tiles <- function(x, window_size) {
  # Make tiles
  tile_borders <- seq(min(x), max(x), window_size)
  # Put each position in its tile
  tiles <- x
  #for (i in seq_along(tiles)) {
  #  tiles[i] <- tail(tile_borders[tiles[i] >= tile_borders], 1)
  #}
  sapply(x, function(x) tail(tile_borders[x >= tile_borders], 1))
}