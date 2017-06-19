#' Collect signal from telomeres for all chromosome arms
#'
#' Pulls out the ChIP signal starting at chromosome ends inward up to a
#' specified distance.\cr
#' \strong{Note:} The fact that some sub-telomeric sequences may be incomplete
#' in the genome means that in some cases we have ChIP-seq data mapping all the
#' way to the very end of the available sequence. This, together with the read
#' extension performed by MACS, leads to some cases where the last positions in
#' the signal data for the right chromosome arms are higher than the annotated
#' length of the respective chromosome. As a result, some negative position
#' values may appear in the output of this function (typically not off by more
#' than 150 bp).\cr
#' \cr \cr
#' @param signal_data Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param length_to_collect Integer specifying the length (in bp) of the region
#' to collect signal for, starting form the telomeres. Defaults to 100000 (i.e.
#' 100 kb).
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' No default.
#' @return An R data frame containing 32 rows (one for each chromosome arm) and
#' a number of columns dependent on the specified \code{length_to_collect}. The
#' following columns are included:
#' \enumerate{
#'   \item \code{chr} The chromosome number
#'   \item \code{arm} The chromosome arm; one of "L" and "R"
#'   \item \code{size-cat} The chromosome size category; one of "small" and
#'   "large"
#'   \item \code{t1:tn} Signal columns for each position; number depends on the
#'   the specified collection length (\code{n = length_to_collect})
#' }
#' @examples
#' \dontrun{
#' signal_from_telomeres2(WT, genome = "SK1Yue")
#' 
#' signal_from_telomeres2(WT, length_to_collect = 50000, genome = "sacCer3")
#' }
#' @export


signal_from_telomeres2 <- function(signal_data, length_to_collect=100000,
                                   genome) {
  t0  <- proc.time()[3]
  
  # IO checks
  check_package("GenomicRanges")
  check_package("EnrichedHeatmap")
  
  if (!is(signal_data, "GRanges")) {
    stop('"signal_data" must be a GRanges object.', call. = FALSE)
  }
  
  if (!is(length_to_collect, "numeric")) {
    stop('"length_to_collect" must be numeric.', call. = FALSE)
  }
  
  if (missing(genome)) stop('"genome" is a required argument.\n', call. = FALSE)
  
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using 'data-raw/data_internal.R'; stored in 'R/sysdata.rda'
  if (genome == 'SK1Yue') {
    coord_table <- SK1Yuecen
  } else if (genome == 'sacCer3') {
    coord_table <- sacCer3cen
  } else if (genome == 'SK1') {
    coord_table <- SK1cen
  } else stop('"genome" argument must be one of "SK1Yue", "sacCer3" or "SK1".')
  
  message('Making GRanges object of subtelomeric regions...')
  # Make sure it is integer
  length_to_collect <- floor(length_to_collect)
  
  # Left arms
  left_arm <- coord_table
  
  GenomicRanges::start(left_arm) <- 1
  GenomicRanges::end(left_arm) <- length_to_collect
  GenomicRanges::strand(left_arm) <- '+'

  # Right arms
  right_arm <- coord_table
  ends <- GenomeInfoDb::seqlengths(right_arm)
  GenomicRanges::end(right_arm) <- ends
  GenomicRanges::start(right_arm) <- (ends - length_to_collect)
  GenomicRanges::strand(right_arm) <- '-'
  
  # Concatenate arms
  telomeres <- c(left_arm, right_arm)
  
  # Compute signal at each gene using package EnrichedHeatmap
  message('Collecting signal...')
  mat <- EnrichedHeatmap::normalizeToMatrix(signal_data, telomeres,
                                            value_column="score",
                                            mean_mode="absolute",
                                            extend=0, k=length_to_collect,
                                            empty_value=NA, smooth=FALSE,
                                            target_ratio=1)
  
  # Prepare final data frame
  df <- as.data.frame(mat)
  size_cat <- ifelse(as.character(telomeres@seqnames) %in% c('chrI',
                                                             'chrIII',
                                                             'chrVI'),
                     'small', 'large')
  
  df <- cbind(chr=as.character(telomeres@seqnames),
              arm=c(rep('L', 16), rep('R', 16)),
              size_cat=size_cat,
              df)
  
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
  return(df)
}