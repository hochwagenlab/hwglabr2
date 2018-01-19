#' Compute spike-in normalization factor
#'
#' Computes spike-in normalization factor between two spiked-in samples. Inputs
#' text files containing counts of aligned reads per chromosome of a hybrid
#' SK1:S288C genome.
#' @param ref_chip_counts Path to reference ChIP samples' read counts file.
#' No default.
#' @param ref_input_counts Path to reference input samples' read counts file.
#' No default.
#' @param test_chip_counts Path to test ChIP samples' read counts file.
#' No default.
#' @param test_input_counts Path to test input samples' read counts file.
#' No default.
#' @return Numeric normalization factor.
#' @examples
#' \dontrun{
#' spikein_normalization_factor(ref_chip_counts='Counts_AH119_chip.txt',
#'                              ref_input_counts='Counts_AH119_input.txt',
#'                              test_chip_counts='Counts_AH8104_chip.txt',
#'                              test_input_counts='Counts_AH8104_input.txt')
#' }
#' @export

spikein_normalization_factor <- function(ref_chip_counts, ref_input_counts,
                                         test_chip_counts, test_input_counts){
  # Put paths in list
  files <- list(ref_chip=ref_chip_counts, ref_input=ref_input_counts,
                test_chip=test_chip_counts, test_input=test_input_counts)
  
  # IO checks
  if (all(unlist(lapply(files, function(x) is(x, 'character'))))) {
      lapply(files, check_path)
  } else stop("arguments must be working paths to input files")
  
  # Print files to read to console
  message('>>> Read alignment count files:')
  for (file in files) {
      message('   ', basename(file))
  }
  
  # Add space to console output
  message()
  # Read files into tibble in list
  tables <- sapply(files, FUN=readr::read_tsv, col_names=F,
                   simplify=FALSE, USE.NAMES=TRUE)
  
  # Add space to console output
  message()
  # Get read counts per chromosome
  message('>>> Count reads per genome:')
  counts <- sapply(tables, FUN=sum_per_genome, simplify=FALSE, USE.NAMES=TRUE)
  
  # Compute normalization factor
  result <- normalization_factor(ctrl_input=counts$ref_input,
                                 ctrl_chip=counts$ref_chip,
                                 test_input=counts$test_input,
                                 test_chip=counts$test_chip)
  
  message('---')
  message('Done!')
  
  return(result)
}

# Helper functions
sum_per_genome <- function(df) {
    # Compute sum of reads aligned to each genome
    S288C <- sum(
        df[apply(df, 1, function(x) stringr::str_detect(x[1],'_S288C')), 2])
    SK1 <- sum(
        df[apply(df, 1, function(x) stringr::str_detect(x[1], '_SK1')), 2])
    
    # Print result to console
    message('  S288C: ', formatC(S288C, big.mark=",",
                                 drop0trailing=TRUE, format="f"))
    message('  SK1: ', formatC(SK1, big.mark=",",
                               drop0trailing=TRUE, format="f"))
    message('      ', round(S288C * 100 / (SK1 + S288C), 1), '% spike-in reads')
    
    # Return result as named vector
    c('S288C'=S288C, 'SK1'=SK1)
}


normalization_factor <- function(ctrl_input, ctrl_chip,
                                 test_input, test_chip) {
    # ChIP vs Input normalization
    ctrl_SK1 <- ctrl_chip['SK1'] / ctrl_input['SK1']
    ctrl_S288C <- ctrl_chip['S288C'] / ctrl_input['S288C']
    
    test_SK1 <- test_chip['SK1'] / test_input['SK1']
    test_S288C <- test_chip['S288C'] / test_input['S288C']
    
    # SK1 vs spike-in normalization
    ctrl <- ctrl_SK1 / ctrl_S288C
    test <- test_SK1 / test_S288C
    
    # Compute and return normalization factor
    test[[1]] / ctrl[[1]]
}
