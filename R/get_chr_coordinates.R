#' Get chromosome coordinates
#'
#' Returns chromosome length and centromere coordinates for specified genome.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' Default is \code{sacCer3}.
#' @param as_df Logical specifying whether the output should be returned as a
#' \code{data frame}. If \code{FALSE}, output is a \code{GRanges} object.
#' Defaults to \code{FALSE}.
#' @return Each chromosome length and centromere start, end and midpoint for each chromosome.
#' @examples
#' \dontrun{
#' get_chr_coordinates()
#' 
#' get_chr_coordinates(genome='SK1', as_df=FALSE)
#' }
#' @export

get_chr_coordinates <- function(genome='sacCer3', as_df=FALSE){
  # IO checks
  check_package("GenomicRanges")
  if (!is(genome, "character")) stop('"genome" must be a character object.')
  
  if (genome == 'sacCer3') {
    coord_table <- sacCer3cen
  } else if (genome == 'SK1') {
    coord_table <- SK1cen
  } else stop('"genome" argument must be either "sacCer3" or "SK1".')
  
  if (as_df) {
    check_package("GenomicRanges")
    coord_table <- cbind(GenomicRanges::as.data.frame(coord_table),
                   chr_len=GenomeInfoDb::seqlengths(coord_table))
  }
    
  return(coord_table)
}