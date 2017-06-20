#' Get Red1 ChIP-seq summits in WT
#'
#' Returns all Red1 ChIP-seq summits called in the average of four WT replicates
#' for which -log10 of the peak's p-value is equal or greater than 20.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#' }
#' Default is \code{"SK1Yue"}.
#' @param as_df Logical specifying whether the output should be returned as a
#' \code{data frame}. If \code{FALSE}, output is a \code{GRanges} object.
#' Defaults to \code{FALSE}.
#' @return Either a \code{Granges} or a {data.frame} object containing Red1
#' summits.
#' @examples
#' \dontrun{
#' get_Red1_summits()
#' 
#' get_Red1_summits(genome='sacCer3', as_df=TRUE)
#' }
#' @export

get_Red1_summits <- function(genome='SK1Yue', as_df=FALSE){
  # IO checks
  check_package("GenomicRanges")
  if (!is(genome, "character")) stop('"genome" must be a character object.')
  
  if (genome == 'SK1Yue') {
    Red1_summits <- SK1Yue_Red1_summits
  } else if (genome == 'sacCer3') {
    Red1_summits <- sacCer3_Red1_summits
  } else if (genome == 'SK1') {
    stop('Data not included for "genome=SK1".\n',
         'Please use either "SK1Yue" or "sacCer3"', call. = FALSE)
  } else stop('"genome" argument must be one of "SK1Yue", "sacCer3" or "SK1".')
  
  if (as_df) {
    Red1_summits <- data.frame(chr=GenomicRanges::seqnames(Red1_summits),
                               start=GenomicRanges::start(Red1_summits),
                               end=GenomicRanges::end(Red1_summits),
                               name=Red1_summits$name,
                               neglog10pvalue=Red1_summits$score)
  }
    
  return(Red1_summits)
}