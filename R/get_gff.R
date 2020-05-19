#' Get GFF
#'
#' Returns GFF table for specified genome.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"S288CYue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' No default.
#' @return GFF table as \code{GRanges} object.
#' @examples
#' \dontrun{
#' get_gff(genome='SK1Yue')
#'
#' get_gff(genome='sacCer3')
#' }
#' @export

get_gff <- function(genome){
  # IO checks
  check_package("GenomicRanges")
  if (!is(genome, "character")) stop('"genome" must be a character object.')

  if (genome == 'SK1Yue') {
      gff <- SK1Yue_gff
  } else if (genome == 'sacCer3') {
      gff <- sacCer3_gff
  } else if (genome == 'SK1') {
      gff <- SK1_gff
  } else if (genome == 'S288CYue') {
      gff <- S288CYue_gff
  } else stop('"genome" argument must be one of "SK1Yue", "S288CYue","sacCer3" or "SK1".')

  return(gff)
}
