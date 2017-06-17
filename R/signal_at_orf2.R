#' Signal at all ORFs genome-wide normalized to constant length (meta ORF)
#'
#' Pulls out the ChIP signal over all ORFs in the genome normalized to a
#' constant length. Collects the signal over each ORF plus both flanking regions
#' (1/2 the length of the ORF on each side) and scales them all to the same
#' value (1000). This means that for two example genes with lengths of 500 bp
#' and 2 kb, flanking regions of 250 bp and 1 kb will be collected up and
#' downstream, respectively. The whole region is then rescaled to a length of
#' 1000, corresponding to a gene length of 500 plus 250 for each flanking
#' region.\cr\cr\cr
#' \strong{Note:} Our ChIP-seq data always contains gaps with missing data. The
#' affected genes will contain "NA" values in the output. As a warning, the
#' number of affected genes is printed to the console.\cr
#' @param signal_data Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param gff Either a path to a gff file or loaded gff data as a \code{GRanges}
#' object. No default.
#' @param write_to_file Logical indicating whether output should be written to a
#' .txt file. The file will be saved at location provided in \code{file_name}
#' argument or at the current working directory if only a file name is provided).
#' If \code{write_to_file = FALSE} the function returns output R object.
#' Defaults to \code{FALSE}.
#' @param file_name Character string indicating file name when writing output to
#' a .txt file. Must be provided if \code{write_to_file = TRUE}. No default.
#' @return An \code{EnrichedHeatmap} matrix with 1000 columns containing signal
#' at each meta position. Columns correspond to the following meta positions:
#' \enumerate{
#'   \item \code{t1 to t250} Upstream flank (1/2 of gene length)
#'   \item \code{t251 to t750} ORF
#'   \item \code{t751 to t1000} Downstream flank (1/2 of gene length)
#' }
#' @examples
#' \dontrun{
#' signal_at_orf2(WT, gff = gff)
#' 
#' signal_at_orf2(WT, gff = "S288C_annotation_modified.gff",
#'                write_to_file = TRUE, file_name = "WT_S288C_metaORF.txt")
#' }
#' @export

signal_at_orf2 <- function(signal_data, gff, write_to_file=FALSE, file_name) {
  t0  <- proc.time()[3]
  
  # IO checks
  check_package("GenomicRanges")
  check_package("EnrichedHeatmap")
  
  if (!is(signal_data, "GRanges")) {
    stop('"signal_data" must be a GRanges object.')
  }
  
  if (missing(gff)) stop('No gff data provided.\n',
                         '"gff" must be gff data as a "GRanges" object\n',
                         'or the path to a gff file.', call. = FALSE)
  
  if (!is(gff, "GRanges")) {
    if (is(gff, "character")) {
      check_path(gff)
      check_package("rtracklayer")
      message('Loading gff file...')
      gff <- rtracklayer::import.gff(gff)
    } else {
      stop('"gff" must be either a GRanges object or a path to a gff file.')
    }
  }
  
  if (write_to_file & missing(file_name)) {
    stop('No file name provided.\n',
         '"file_name" must be be provided when "write_to_file = TRUE"',
         call. = FALSE)
  }
  
  # Drop 'chrMito' and '2-micron' if present in gff (absent from ChIP-seq data)
  gff <- GenomeInfoDb::dropSeqlevels(gff, c('chrMito', '2-micron'))
  
  # Check seqnames (must match between input data and gff)
  # this does not catch the problem of mixing SK1Yue and S288C!
  if (check_chr_names(signal_data) != check_chr_names(gff)) {
    stop("The reference genomes in the data and the gff do not seem to match.",
         call. = FALSE)
  } else if (check_chr_names(signal_data) == check_chr_names(gff)) {
    message('Chromosomes numbered using ', check_chr_names(signal_data))
  } else stop('Did not recognize reference genome.\n',
              'Please ensure chromosome numbers are in the expected format:\n',
              'e.g. "chrI" or "chr01".')
  
  message('Types of features in the gff data:')
  for(i in 1:length(unique(gff$type))) {
    message(paste0('   ', unique(gff$type)[i]))
  }
  
  # Add 1/2-of-gene-length flanks to genes
  # (can ignore strand since it is the same on both ends)
  message('Preparing gene flanks...')
  flank <- floor(GenomicRanges::width(gff) / 2)
  GenomicRanges::start(gff) <- GenomicRanges::start(gff) - flank
  GenomicRanges::end(gff) <- GenomicRanges::end(gff) + flank
  # Fix overflowing flanks (ensure they don't go over chromosome bounds)
  gff <- GenomicRanges::trim(gff)
  
  # Compute signal at each gene using package EnrichedHeatmap
  message('Computing signal at each normalized ORF...')
  mat <- EnrichedHeatmap::normalizeToMatrix(signal_data, gff,
                                            value_column="score",
                                            mean_mode="absolute",
                                            extend=0, k=1000, empty_value=NA,
                                            smooth=FALSE, target_ratio=1)
  
  # Check NAs
  n_genes <- nrow(mat)
  n_complete_genes <- nrow(mat[complete.cases(mat), ])
  n_genes_with_NAs <- n_genes - n_complete_genes
  message('---')
  message('Genes containing NA values: ', n_genes_with_NAs,
          ' of a total of ', n_genes, ' (',
          round(n_genes_with_NAs / n_genes * 100, 1), '%)')

  if (write_to_file){
    message(paste0('Writing data to file: ', file_name))
    write.table(mat, file_name, sep="\t", quote=FALSE, row.names=FALSE)
    
    message('---')
    message('Completed in ', hwglabr2::elapsed_time(t0))
  } else {
    # Adjust attributes of EnrichedHeatmap matrix:
    attr(mat, "upstream_index") <- 1:250
    attr(mat, "target_index") <- 251:750
    attr(mat, "downstream_index") <- 751:1000
    attr(mat, "extend") <- c('', '') # Leave x axis 1 and 1000 annotations blank
    attr(mat, "signal_name") <- deparse(substitute(signal_data))
    attr(mat, "target_name") <- 'ORFs'
    
    message('---')
    message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
    return(mat)
    }
}
