#' Compute correlation between signal track data for two samples
#'
#' Returns correlation between ChIP-seq signal of two samples. Useful to check
#' replicates or compare different samples.
#' 
#' @param signal_data_A Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param signal_data_B Second samples' signal track data in the same format as
#' \code{signal_data_A} object. No default.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' No default.
#' @param method Character object specifying which correlation coefficient is to
#' be computed. Will be the input to the \code{method} argument of function
#' \code{cor} of the \code{stats} package. Accepts one of (can be abbreviated):
#' \enumerate{
#'   \item \code{"pearson"}
#'   \item \code{"kendall"}
#'   \item \code{"spearman"}
#' }
#' Defaults to \code{"pearson"}.
#' @return A floating point value corresponding to the signal correlation.
#' @examples
#' \dontrun{
#' signal_track_correlation(WT, dot1, genome = "SK1Yue")
#' 
#' signal_track_correlation(WT_A, WT_B, genome = "sacCer3", method="spearman")
#' }
#' @export

signal_track_correlation <- function(signal_data_A, signal_data_B, genome,
                                     method="pearson") {
  t0  <- proc.time()[3]
  
  # IO checks
  check_package("GenomicRanges")
  
  if (!is(signal_data_A, "GRanges") | !is(signal_data_B, "GRanges")) {
    stop('"signal_data_A" and "signal_data_B" must be GRanges objects.',
         call. = FALSE)
  }
  
  if (missing(genome)) stop('"genome" is a required argument.\n', call. = FALSE)
  
  # Get "seqlengths"
  message('Get genome coordinates...')
  genome_info <- hwglabr2::get_chr_coordinates(genome=genome)
  
  # Sort sequences and levels to make sure they match
  genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
  grA <- sort(GenomeInfoDb::sortSeqlevels(signal_data_A))
  grB <- sort(GenomeInfoDb::sortSeqlevels(signal_data_B))
  
  # Add info to signal object
  suppressWarnings(
    GenomeInfoDb::seqlengths(grA) <- GenomeInfoDb::seqlengths(genome_info)
  )
  suppressWarnings(
    GenomeInfoDb::seqlengths(grB) <- GenomeInfoDb::seqlengths(genome_info)
  )
  
  # Compute 1-bp tiling windows
  message('Compute 1-bp tiling windows...')
  bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(grA),
                                    tilewidth=1,
                                    cut.last.tile.in.chrom=TRUE)
  
  # Get signal as "RleList"; the signal is stored in the "score" metadata column
  message('Compute signal per bp (or tile)...')
  grA <- GenomicRanges::coverage(grA, weight="score")
  grB <- GenomicRanges::coverage(grB, weight="score")
  
  # Get signal per tile
  A_bins <- GenomicRanges::binnedAverage(bins, grA, "binned_score")
  B_bins <- GenomicRanges::binnedAverage(bins, grB, "binned_score")
  
  # Replace zeros by R's representation of missing data
  message('Replace zeros by "NA"s...')
  A_bins[A_bins$binned_score == 0]$binned_score <- NA
  B_bins[B_bins$binned_score == 0]$binned_score <- NA
  
  message('Compute ', method, ' correlation between samples: ',
          deparse(substitute(signal_data_A)), ' and ',
          deparse(substitute(signal_data_B)), '...')
  corr <- cor(x=A_bins$binned_score, y=B_bins$binned_score,
              use="complete.obs", method=method)
  
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
  
  corr
}