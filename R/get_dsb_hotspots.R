#' Get DSB hotspots
#'
#' Returns double-strand break (DSB) hotspots mapped using the Spo11 oligo
#' technique (hotspot location data from 
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/21376234}{Pan \emph{et al.} 2011};
#' hotspot intensities from 
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/24717437}{Thacker \emph{et al.} 
#' 2014}).
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#' }
#' No default.
#' @param as_df Logical specifying whether the output should be returned as a
#' \code{data frame}. If \code{FALSE}, output is a \code{GRanges} object.
#' Defaults to \code{FALSE}.
#' @return Either a \code{Granges} or a {data.frame} object containing DSB
#' hotspot start and end coordinates, as well as intensity.
#' @examples
#' \dontrun{
#' get_dsb_hotspots(genome='SK1Yue')
#' 
#' get_dsb_hotspots(genome='sacCer3', as_df=TRUE)
#' }
#' @export

get_dsb_hotspots <- function(genome, as_df=FALSE){
  # IO checks
  check_package("GenomicRanges")
  if (!is(genome, "character")) stop('"genome" must be a character object.')
  
  if (genome == 'SK1Yue') {
    Spo11_DSBs <- SK1Yue_Spo11_DSBs
  } else if (genome == 'sacCer3') {
    Spo11_DSBs <- sacCer3_Spo11_DSBs
  } else if (genome == 'SK1') {
    stop('Data not included for "genome=SK1".\n',
         'Please use either "SK1Yue" or "sacCer3"', call. = FALSE)
  } else stop('"genome" argument must be one of "SK1Yue" or "sacCer3".')
  
  if (as_df) {
    Spo11_DSBs <- data.frame(chr=GenomicRanges::seqnames(Spo11_DSBs),
                             start=GenomicRanges::start(Spo11_DSBs),
                             end=GenomicRanges::end(Spo11_DSBs),
                             score=Spo11_DSBs$score)
  }
    
  return(Spo11_DSBs)
}