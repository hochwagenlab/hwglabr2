#' Get intergenic regions
#'
#' Returns table of intergenic regions for specified genome.
#' @param genome Character object specifying the reference genome; accepts one
#' of the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' No default.
#' @param as_gr Logical specifying whether the output should be returned as a
#' \code{GRanges} object. If \code{FALSE}, output is a \code{data.frame}. Can
#' only be used if \code{genome="SK1Yue"}, because \code{"SK1"} and
#' \code{"sacCer3"} contain overlapping genes, resulting in negative intervals
#' (not allowed by \code{IRanges}). Defaults to \code{FALSE}.
#' @return R data frame containing coordinates and annotations for all
#' intergenic regions in the genome.
#' @examples
#' \dontrun{
#' get_intergenic_regions(genome='SK1Yue')
#' 
#' get_intergenic_regions(genome='sacCer3')
#' 
#' }
#' @export

get_intergenic_regions <- function(genome, as_gr=FALSE){
  # IO checks
  if (as_gr & genome != 'SK1Yue') {
    stop('"as_gr" can only be used when genome="SK1Yue".', call. = FALSE)
  }
  if (!is(genome, "character")) stop('"genome" must be a character object.')
  
  if (genome == 'SK1Yue') {
    inter_table <- SK1Yue_intergenic
  } else if (genome == 'sacCer3') {
    inter_table <- sacCer3_intergenic
  } else if (genome == 'SK1') {
    inter_table <- SK1_intergenic
  } else stop('"genome" argument must be one of "SK1Yue", "sacCer3" or "SK1".')
  
  if (as_gr) {
    inter_table <- GenomicRanges::GRanges(
      seqnames=inter_table$chr,
      ranges=IRanges::IRanges(start=inter_table$left_coordinate,
                              end=inter_table$right_coordinate),
      strand="+", type=inter_table$type,
      left_gene=inter_table$left_gene,
      left_gene_strand=inter_table$left_gene_strand,
      right_gene=inter_table$right_gene,
      right_gene_strand=inter_table$right_gene_strand
    )
  }

  return(inter_table)
}