#' Get chromosome coordinates
#'
#' Returns chromosome length and centromere coordinates for specified genome.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' Default is \code{"SK1Yue"}.
#' @param as_df Logical specifying whether the output should be returned as a
#' \code{data frame}. If \code{FALSE}, output is a \code{GRanges} object.
#' Defaults to \code{FALSE}.
#' @return Centromere start and end, as well as chroosome lengths, for each
#' chromosome.
#' @examples
#' \dontrun{
#' get_chr_coordinates()
#' 
#' get_chr_coordinates(genome='sacCer3', as_df=FALSE)
#' 
#' get_chr_coordinates(genome='SK1', as_df=TRUE)
#' }
#' @export

get_chr_coordinates <- function(genome='SK1Yue', as_df=FALSE){
  # IO checks
  check_package("GenomicRanges")
  if (!is(genome, "character")) stop('"genome" must be a character object.')
  
  if (genome == 'SK1Yue') {
    coord_table <- SK1Yuecen
  } else if (genome == 'sacCer3') {
    coord_table <- sacCer3cen
  } else if (genome == 'SK1') {
    coord_table <- SK1cen
  } else stop('"genome" argument must be one of "SK1Yue", "sacCer3" or "SK1".')
  
  if (as_df) {
    coord_table <- cbind(GenomicRanges::as.data.frame(coord_table),
                   chr_len=GenomeInfoDb::seqlengths(coord_table))
  }
    
  return(coord_table)
}