#' Signal flanking rDNA
#'
#' Pulls out ChIP signal at specified window sizes flanking the rDNA region on
#' chromosome 12.
#' @param signal_data Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param flank_length Integer specifying the length (in bp) of the windows to
#' collect signal for up and downstream of the rDNA. Defaults to 40000 (i.e.
#' 40 kb).
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' No default.
#' @examples
#' \dontrun{
#' signal_flanking_rDNA(WT, genome = 'SK1Yue')
#' 
#' signal_flanking_rDNA(WT, flank_length=50000, genome = 'sacCer3')
#' }
#' @export

signal_flanking_rDNA <- function(signal_data, flank_length=40000, genome) {
    # IO checks
    check_package("GenomicRanges")
    
    if (!is(signal_data, "GRanges")) {
        stop('"signal_data" must be a GRanges object.')
    }
    
    if (missing(genome)) stop('"genome" is a required argument.\n', call. = FALSE)
    
    message('\nCollecting signal...')
    
    ### rDNA coordinates and keep only chrXII
    if (genome == 'SK1Yue') {
        start <- 447012 - flank_length
        end <- 461699 + flank_length
        signal_data <- signal_data[signal_data@seqnames == 'chrXII']
    } else if (genome == 'sacCer3') {
        start <- 451575 - flank_length
        # A stretch present in S288C downstream of rDNA is absent in SK1
        # Adapt coordinates accordingly (until about bp 490'500)
        #end <- 468931 + flank_length
        end <- 490500 + flank_length
        signal_data <- signal_data[signal_data@seqnames == 'chrXII']
    } else if (genome == 'SK1') {
        start <- 433029 - flank_length
        end <- 451212 + flank_length
        signal_data <- signal_data[signal_data@seqnames == 'chr12']
    } else stop('"genome" argument must be one of "SK1Yue", "sacCer3" or "SK1".')
    
    # Convert positions to indices
    start <- tail(signal_data[signal_data@ranges@start < start, 1], 1)
    start <- which(signal_data@ranges@start == start@ranges@start)
    
    end <- head(signal_data[signal_data@ranges@start > end, 1], 1)
    end <- which(signal_data@ranges@start == end@ranges@start)
    
    # Get subset of data corresponding to rDNA +/- flank_length
    rDNA_signal <- signal_data[start:end]
    
    message('Done!')
    return(rDNA_signal)
}
