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

get_intergenic_regions <- function(genome){
  # IO checks
  if (!is(genome, "character")) stop('"genome" must be a character object.')
  
  if (genome == 'SK1Yue') {
    inter_table <- SK1Yue_intergenic
  } else if (genome == 'sacCer3') {
    inter_table <- sacCer3_intergenic
  } else if (genome == 'SK1') {
    inter_table <- SK1_intergenic
  } else stop('"genome" argument must be one of "SK1Yue", "sacCer3" or "SK1".')

  return(inter_table)
}