#' Signal at all ORFs genome-wide normalized to constant length (meta ORF)
#'
#' This function allows you to pull out the ChIP signal over all ORFs in the
#' genome (normalized to a constant length). It collects the signal over each
#' ORF plus both flanking regions (1/2 the length of the ORF on each side) and
#' scales them all to the same value (1000). This means that for two example
#' genes with lengths of 500 bp and 2 kb, flanking regions of 250 bp and 1 kb
#' will be collected up and downstream, respectively. The whole region is then
#' rescaled to a length of 1000, corresponding to a gene length of 500 plus 250
#' for each flanking region.\cr
#' \cr \cr
#' \strong{Note:} Our ChIP-seq data always contains gaps with missing data. The
#' affected genes will contain "NA" values in the output. As a warning, the
#' number of affected genes is printed to the console.\cr
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param gff Either a path to a gff file or loaded gff data as a \code{GRanges}
#' object. No default.
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
#' signal_at_orf2(WT, gff = gff)
#' 
#' signal_at_orf2(WT, gff = "S288C_annotation_modified.gff",
#'                write_to_file = TRUE)
#' }
#' @export

signal_at_orf2 <- function(gr, gff, write_to_file=FALSE) {
  t0  <- proc.time()
  
  # IO checks
  check_package("GenomicRanges")
  check_package("EnrichedHeatmap")
  
  if (!is(gr, "GRanges")) stop('"gr" must be a GRanges object.')
  
  if (missing(gff)) stop('No gff data provided.\n',
                         '"gff" must be gff data as a "GRanges" object\n',
                         'or the path to a gff file.', call. = FALSE)
  
  if (!is(gff, "GRanges") && is(gff, "character")) {
    check_path(gff)
    check_package("rtracklayer")
    message('Loading gff file...')
    gff <- rtracklayer::import.gff(gff)
  } else stop('"gff" must be either a GRanges object or a path to a gff file.')
  
  # Drop 'chrMito' and '2-micron' if present in gff (absent from ChIP-seq data)
  gff <- gff[!as.character(gff@seqnames) %in% c('chrMito', '2-micron'), ]
  
  # Check reference genome (must match between input data and gff)
  if (check_genome(gr)[1] != check_genome(gff)[1]) {
    stop("The reference genomes in the data and the gff do not seem to match.",
         call. = FALSE)
  } else if (check_genome(gr)[1] == check_genome(gff)[1]) {
    message('Ref. genome: ', paste(check_genome(gr), collapse = " "))
  } else stop('Did not recognize reference genome.\n',
              'Please ensure chromosome numbers are in the expected format:\n',
              'e.g. "chrI" or "chr01".')
  
  message('Types of features in the gff data:')
  for(i in 1:length(unique(gff$type))) {
    message(paste0('   ', unique(gff$type)[i]))
  }
  
  # Add 1/2-of-gene-length flanks to genes
  message('Preparing gene flanks...')
  flank <- floor(GenomicRanges::width(gff) / 2)
  GenomicRanges::start(gff) <- GenomicRanges::start(gff) - flank
  GenomicRanges::end(gff) <- GenomicRanges::end(gff) + flank
  # Fix overflowing flanks?
  GenomicRanges::start(gff)[GenomicRanges::start(gff) < 0] <- 1
  # Add chr lengths to fix the other way?
  
  # Compute signal at each gene using package EnrichedHeatmap
  message('Computing signal at each normalized ORF...')
  mat <- EnrichedHeatmap::normalizeToMatrix(gr, gff, value_column="score",
                                            mean_mode="absolute",
                                            extend=0, k=1000, empty_value=NA,
                                            smooth=FALSE, target_ratio=1)
  
  ### Add gene names to data
  # GFF files are not guaranteed to contain a consistent column name for gene ID
  # Check columns in provided GFF in search for the gene ID column
  if ('group' %in% names(GenomicRanges::mcols(gff))) {
    df <- cbind("gene"=gff$group, as.data.frame(mat))
  } else if ('ID' %in% names(GenomicRanges::mcols(gff))) {
    df <- cbind("gene"=gff$ID, as.data.frame(mat))
  } else {
    message('Could not find gene IDs in the provided gff')
    message('(the output will not include gene IDs).')
    df <- as.data.frame(mat)
  }
  
  # Check NAs
  n_genes <- nrow(df)
  n_complete_genes <- nrow(df[complete.cases(df), ])
  n_genes_with_NAs <- n_genes - n_complete_genes
  message('---')
  message('Genes containing NA values: ', n_genes_with_NAs,
          ' of a total of ', n_genes, ' (',
          round(n_genes_with_NAs / n_genes * 100, 1), '%)')
  
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0))
  
  if (write_to_file){
    file_name <- paste0(deparse(substitute(gr)), "_",
                        check_genome(gr)[1], "_metaORF.txt")
    message(paste0('Writing data to file: ', file_name))
    write.table(df, file_name, sep="\t", quote=FALSE, row.names=FALSE)
    
    message('Done!')
  } else return(df)
}
