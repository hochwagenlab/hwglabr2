#' Average signal genome-wide and on each chromosome
#'
#' Computes average signal (in the \code{score} \code{GRanges} metadata column)
#' genome-wide and on each chromosome (each individual sequence determined by
#' \code{seqnames}).
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param remove_cen Logical indicating whether to remove regions around
#' centromeres. Defaults to \code{FALSE}.
#' @param genome Character string indicating reference genome used to align
#' the data. Must be provided when \code{remove_cen = TRUE}, in order to load
#' appropriate centromere data. Accepts one of the following strings:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#'   \item \code{"SK1"}
#' }
#' No default.
#' @param cen_region_length Integer indicating the length (in bp) of the region
#' to remove (centered on the centromere of each chromosome).
#' Defaults to 50'000 bp.
#' @param mean_norm Logical indicating whether to normalize by average
#' genome-wide score. This adds two columns to the dataframe (first element of
#' the output list): mean ratio (\code{avrg_signal_mean_ratio}) and
#' mean-subtracted (\code{avrg_signal_mean_subtracted}).
#' Defaults to \code{FALSE}.
#' @param order_chrs Logical indicating whether to order rows by chromosome
#' length (in increasing order). Defaults to \code{FALSE}.
#' @return List with two elements:
#' \enumerate{
#'   \item \code{seq_avrg} Average \code{score} on each sequence (dataframe)
#'   \item \code{genome_avrg} Average \code{score} genome-wide (named vector)
#' }
#' @examples
#' \dontrun{
#' average_chr_signal(GRanges_object)
#' 
#' average_chr_signal(anti_Rec8, remove_cen=TRUE, mean_norm=TRUE, order_chrs=TRUE)
#' 
#' average_chr_signal(anti_Rec8, remove_cen=TRUE, genome = 'SK1Yue',
#'                    cen_region_length=40000, mean_norm=TRUE, order_chrs=TRUE)
#' }
#' @export

average_chr_signal <- function(gr, remove_cen=FALSE, genome,
                               cen_region_length=50000,
                               mean_norm = FALSE, order_chrs=FALSE){
  # IO checks
  check_package("GenomicRanges")
  if (!is(gr, "GRanges")) stop('input must be a GRanges object.')
  if (!"score" %in% names(GenomicRanges::mcols(gr))) {
    stop(deparse(substitute(gr)), ' does not have a "score" metadata column.')
  }
  
  if (remove_cen) {
    if (missing(genome)) {
      stop('No reference genome provided.\n',
           '"genome" is required when "remove_cen = TRUE"',
           call. = FALSE)
    }
    
    if (!is(genome, "character")) stop('"genome" must be a character object.')
    
    # Make sure provided genome does not clash with data
    if (check_chr_names(gr) == 'roman numerals') {
      if (!genome %in% c('SK1Yue', 'S288C')) {
        stop('"genome" must be either "SK1Yue" or "S288C"', call. = FALSE)
      }
    } else {
      if (genome != 'SK1') stop('"genome" must be "SK1"', call. = FALSE)
    }
    
    message('Removing ', floor(cen_region_length / 1000),
            '-Kb regions around centromeres...')

    # Load centromere data
    if (genome == 'SK1Yue') {
      cen <- SK1Yuecen
    } else if (genome == 'sacCer3') {
      cen <- sacCer3cen
    } else if (genome == 'SK1') {
      cen <- SK1cen
    }
    
    # Add/subtract half of region length (centered on centromere midpoint)
    half_length <- floor(cen_region_length / 2)
    offset <- floor(GenomicRanges::width(cen) / 2)
    
    GenomicRanges::start(cen) <- (GenomicRanges::start(cen) + offset
                                  - half_length)
    GenomicRanges::end(cen) <- (GenomicRanges::end(cen) - offset
                                + half_length)
    # Remove centromere regions
    gr <- gr[!GenomicRanges::overlapsAny(gr, cen)]
  }
  
  message('Computing average signal...')
  avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x))
                       / sum(GenomicRanges::width(x)))
  genome_avrg <- avrg(gr)
  seq_avrg <- sapply(GenomicRanges::split(gr, GenomicRanges::seqnames(gr)),
                     avrg)
  
  # Convert to dataframe
  seq_avrg <- data.frame(chr=names(seq_avrg), avrg_signal=seq_avrg,
                         row.names=NULL, stringsAsFactors = F)
  
  if (mean_norm) {
    seq_avrg$avrg_signal_mean_subtracted <- seq_avrg$avrg_signal - genome_avrg
    seq_avrg$avrg_signal_mean_ratio <- seq_avrg$avrg_signal / genome_avrg
  }
  
  if (order_chrs) seq_avrg <- order_chromosomes(seq_avrg, chr_column='chr',
                                                decreasing=FALSE)
  
  message('Done!')
  return(list("seq_avrg"=seq_avrg, "genome_avrg"=genome_avrg))
}


order_chromosomes <- function(df, chr_column, decreasing=FALSE) {
  # Make sure input is a data frame with 16 rows
  if (!is.data.frame(df) | nrow(df) != 16) {
    stop("Wrong input data - not an R data frame wit 16 rows.", call. = FALSE)
  }
  
  # Check reference genome and get order
  chrom_roman <- c('chrI', 'chrVI', 'chrIII', 'chrIX', 'chrVIII', 'chrV',
                   'chrXI', 'chrX', 'chrXIV', 'chrII', 'chrXIII', 'chrXVI',
                   'chrXII', 'chrVII', 'chrXV', 'chrIV')
  chrom_arabic <- c('chr01', 'chr06', 'chr03', 'chr09', 'chr08', 'chr05',
                    'chr11', 'chr10', 'chr14', 'chr02', 'chr13', 'chr16',
                    'chr12', 'chr07', 'chr15', 'chr04')
  
  check_roman <- any(grep('chr[XVI]', df[, chr_column]))
  check_arabic <- any(grep('chr[0-9]', df[, chr_column]))
  
  if (check_roman && !check_arabic) {
    chrom <- chrom_roman    
  } else if (check_arabic && !check_roman) {
    chrom <- chrom_arabic
  } else stop('Did not recognize chromosome numbering system.\n',
              'Please make sure chromosomes are numbered according to one of',
              ' the following formats:\n',
              'e.g. "chrI" or "chr01"', call. = FALSE)
  
  if (decreasing) chrom <- rev(chrom)
  
  return(df[match(chrom, df[[chr_column]]), ])
}
