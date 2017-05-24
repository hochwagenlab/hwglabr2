#' Signal around provided genomic positions
#'
#' Pulls out signal in intervals of the defined length centered at provided
#' positions.\cr
#' \cr \cr
#' \strong{Note:} Our ChIP-seq data always contains gaps with missing data. The
#' affected intervals will contain "NA" values in the output. As a warning, the
#' number of affected genes is printed to the console.\cr
#' @param signal_data Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param positions Position data as a \code{GRanges} object. Will also accept
#' an interval \code{GRanges} object of constant width (extended by specified
#' \code{up_ext} and \code{down_ext}). No default.
#' @param position_names String indicating the column containing the position
#' names: either the \code{GRanges} 'seqnames' or a metadata column. Will be
#' added to the output to identify the different positions. No default.
#' @param up_ext Integer indicating number of bp to collect upstream of start
#' position. Defaults to \code{1000}.
#' @param down_ext Integer indicating number of bp to collect downstream of end
#' position. Deafaults to \code{1000}.
#' @param write_to_file Logical indicating whether output should be written to a
#' .txt file (in current working directory). If \code{to_file = FALSE} the
#' function returns output R object. Defaults to \code{FALSE}.
#' @return A data frame with 1001 columns:
#' \enumerate{
#'   \item \code{gene} First column with the gene IDs (if found in input gff)
#'   \item \code{t1 to t1000} 1000 columns with signal at each meta position.
#' }
#' Columns t1 to t1000 correspond to the following meta positions:
#' \enumerate{
#'   \item \code{t1 to t250} Upstream flank (1/2 of gene length)
#'   \item \code{t251 to t750} ORF
#'   \item \code{t751 to t1000} Downstream flank (1/2 of gene length)
#' }
#' @examples
#' \dontrun{
#' signal_at_position(WT, Red1_summits, 'name')
#' 
#' signal_at_position(WT, centromeres, 'seqnames',
#'                    up_ext = 5000, down_ext = 5000,
#'                    write_to_file = TRUE)
#' }
#' @export

signal_at_position <- function(signal_data, positions, position_names,
                               up_ext=1000, down_ext=1000,
                               write_to_file=FALSE) {
  t0  <- proc.time()
  
  # IO checks
  if (missing(signal_data)) stop('No signal data provided: ',
                                 '"signal_data" argument is required',
                                 call. = FALSE)
  if (missing(positions)) stop('No position data provided: ',
                               '"positions" argument is required',
                               call. = FALSE)
  if (missing(position_names)) stop('No position data provided: ',
                                     '"position_names" argument is required',
                                     call. = FALSE)
  
  check_package("GenomicRanges")
  check_package("EnrichedHeatmap")
  
  if (!is(signal_data, "GRanges")) stop('"signal_data" must be a GRanges object.')
  if (!is(positions, "GRanges")) stop('"positions" must be a GRanges object')
  if (!is(position_names, "character")) stop('"position_names" must be an ',
                                              'object of class "character"')
  if (!position_names %in% c('seqnames',
                              names(GenomicRanges::mcols(positions)))) {
    stop('"position_names" must be either "seqnames" or a GRanges metadata ',
         'column name')
  }
  
  # Check reference genome (must match between signal and positions data)
  if (check_genome(signal_data)[1] != check_genome(positions)[1]) {
    stop("The reference genomes in the signal and positions data do not match.",
         call. = FALSE)
  } else if (check_genome(signal_data)[1] == check_genome(positions)[1]) {
    message('Ref. genome: ', paste(check_genome(signal_data), collapse = " "))
  } else stop('Did not recognize reference genome.\n',
              'Please ensure chromosome numbers are in the expected format:\n',
              'e.g. "chrI" or "chr01".')
  
  # Number of positions to analyze
  message('Number of genomic positions:')
  all_pos <- 0
  for (chr in unique(positions@seqnames)) {
    number_pos <- length(positions[positions@seqnames == chr])
    message('     ', chr, ':  ', number_pos)
    all_pos <- all_pos + number_pos
  }
  message('     ...')
  message('     Total of ', all_pos)
  
  # Drop 'chrMito' and '2-micron' if present in positions
  # (absent from ChIP-seq data)
  if (length(GenomicRanges::seqinfo(positions)) == 18) {
    positions <- GenomeInfoDb::dropSeqlevels(positions, c('chrMito', '2-micron'))
  }
  
  # Add chromosome seqlengths (if not present)
  if (any(is.na(seqlengths(positions)))) {
    genome <- ifelse(check_genome(gff)[1] == 'S288c', 'sacCer3', 'SK1')
    chr_lengths <- GenomeInfoDb::seqlengths(hwglabr2::get_chr_coordinates(genome=genome))
    GenomeInfoDb::seqlengths(positions) <- chr_lengths
  }
  
  # Add extensions to intervals
  message('Preparing intervals...')
  up_ext <- up_ext + 1
  #GenomicRanges::start(positions) <- GenomicRanges::start(positions) - up_ext
  #GenomicRanges::end(positions) <- GenomicRanges::end(positions) + down_ext
  # Accounting for strand!
  positions <- suppressWarnings(GenomicRanges::promoters(positions,
                                                         upstream = up_ext,
                                                         downstream = down_ext))
  # Fix overflowing flanks
  # (ensure flanks do not go over chromosome bounds)
  positions <- GenomicRanges::trim(positions)
  
  # Compute signal at each interval using package EnrichedHeatmap
  message('Collecting signal...')
  output <- EnrichedHeatmap::normalizeToMatrix(signal_data, positions,
                                               value_column="score",
                                               mean_mode="absolute", extend=0,
                                               k=up_ext + down_ext,
                                               empty_value=NA,
                                               smooth=FALSE, target_ratio=1)
  
  # Add position names to data as data frame
  if (position_names == 'seqnames') {
    name <- positions@'seqnames'
  } else name <- GenomicRanges::mcols(positions)[position_names]
  
  output <- data.frame('name'=name, output)
  
  # Check NAs
  n_int <- nrow(output)
  n_complete_int <- nrow(output[complete.cases(output), ])
  n_int_with_NAs <- n_int - n_complete_int
  message('---')
  message('Intervals containing NA values: ', n_int_with_NAs,
          ' of a total of ', n_int, ' (',
          round(n_int_with_NAs / n_int * 100, 1), '%)')
  
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0))
  
  if (write_to_file){
    file_name <- paste0(deparse(substitute(signal_data)), '_',
                        check_genome(signal_data)[1], '_signal_at_',
                        deparse(substitute(positions)), '.txt')
    message(paste0('Writing data to file: ', file_name))
    write.table(output, file_name, sep="\t", quote=FALSE, row.names=FALSE)
    
    message('Done!')
  } else return(output)
}
