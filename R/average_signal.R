#' Average signal genome-wide and on each chromosome
#'
#' Computes average signal (in the \code{score} \code{GRanges} metadata column)
#' genome-wide and on each chromosome (each individual sequence determined by
#' \code{seqnames}).
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @return List of two named numeric vectors with average \code{score} on:
#' \enumerate{
#'   \item \code{seq_avrg} Each individual sequence
#'   \item \code{genome_avrg} Whole genome
#' }
#' @examples
#' \dontrun{
#' average_signal(GRanges_objct)
#' 
#' average_signal(WT)
#' }
#' @export

average_signal <- function(gr){
  # IO checks
  check_package("GenomicRanges")
  if (!is(gr, "GRanges")) stop('input must be a GRanges object.')
  if (!"score" %in% names(GenomicRanges::mcols(gr))) {
    stop(deparse(substitute(gr)), ' does not have a "score" metadata column.')
  }
  
  message('Computing average signal...')
  avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x))
                      / sum(GenomicRanges::width(x)))
  genome_avrg <- avrg(gr)
  seq_avrg <- sapply(GenomicRanges::split(gr,
                                          GenomicRanges::seqnames(gr)), avrg)
  
  message('Done!')
  return(list("seq_avrg"=seq_avrg, "genome_avrg"=genome_avrg))
}