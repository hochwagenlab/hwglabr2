#' Apply function to bins of a numerical variable defined along a genome
#'
#' This is a generalization of \code{GenomicRanges::binnedAverage} to allow the
#' use of any applicable function instead of simply the mean. The interface is
#' the same with the added argument for the function to use.
#' @param bins \code{GRanges} object representing the genomic bins. Typically
#' obtained by calling \code{tileGenome} with
#' \code{cut.last.tile.in.chrom=TRUE}. No default.
#' @param numvar Named RleList object representing a numerical variable defined
#' along the genome covered by bins (which is the genome described by
#' \code{seqinfo(bins)}). No default.
#' @param varname Character string indicating the name of the genomic variable.
#' No default.
#' @param fun R function to apply to the numerical variable in each bin.
#' No default.
#' @param ... optional arguments to \code{fun}.
#' @examples
#' \dontrun{
#' # Generate 100-bp genome tiles
#' bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(gr),
#'                                   tilewidth=100, cut.last.tile.in.chrom=TRUE)
#' 
#' # Get signal as "RleList"; stored in the "score" metadata column
#' score <- GenomicRanges::coverage(gr, weight="score")
#' 
#' # Now use the function to get the maximum score in each tile
#' max_per_bin <- binned_function(bins=bins, numvar=score,
#'                                varname="binned_score", fun=max)
#' }
#' @export

binned_function <- function(bins, numvar, varname, fun, ...) {
  t0  <- proc.time()[3]
  
  # IO checks
  check_package("GenomicRanges")
  if (!is(bins, "GRanges")) stop('"bins" must be a GRanges object.')
  if (!is(numvar, "RleList")) stop('"numvar" must be an RleList object.')
  if (!identical(GenomeInfoDb::seqlevels(bins), names(numvar))) {
    stop('seqlevels must match between "bins" and "numvar".')
  }
  
  bins_per_chrom <- split(IRanges::ranges(bins),
                          GenomeInfoDb::seqnames(bins))
  fun_list <- lapply(names(numvar),
                     function(seqname) {
                       views <- IRanges::Views(numvar[[seqname]],
                                               bins_per_chrom[[seqname]])
                       IRanges::viewApply(views, fun)
                     })
  new_mcol <- unsplit(fun_list, as.factor(GenomeInfoDb::seqnames(bins)))
  
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
  return(new_mcol)
}