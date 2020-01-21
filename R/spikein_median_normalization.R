#' Spike-in normalization function
#'
#' Computes median signal (in the \code{gmedian} \code{GRanges} metadata column)
#' genome-wide.
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_bedGraph}}. No default.
#' @return GRanges object with \code{spikescore} median-subtracted,
#' spike-in normalized \code{score} on each sequence
#' }
#' @examples
#' \dontrun{
#' spikein_median_normalization(GRanges_object)
#' }
#' @export

spikein_median_normalization <- function(gr,spikein_factor){
  # IO checks
  if (!is(gr, "GRanges")) stop('input must be a GRanges object.')
  if (!"score" %in% names(GenomicRanges::mcols(gr))) {
    stop(deparse(substitute(gr)), ' does not have a "score" metadata column.')
  }

  message('Computing genome median...')
  gmedian <- function(x) (median(rep(GenomicRanges::score(gr),
                      GenomicRanges::width(gr))))
  genome_median <- gmedian(gr)
  print(genome_median)

  gr$spikescore <- gr$score - genome_median
  gr$spikescore <- gr$spikescore * spikein_factor

  message('Done!')
  return(gr)
}
