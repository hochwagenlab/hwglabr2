# Utility functions used by hwglabr2
#
# Set of utility functions internal to the package that perform auxiliary tasks
# used repeatedly by `hwglabr2` functions.

check_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package == "EnrichedHeatmap") {
      stop('Requires R package "EnrichedHeatmap". Please install it.',
           'Note: please install the latest version ',
           '(from GitHub, NOT Bioconductor):\n',
           'devtools::install_github("jokergoo/EnrichedHeatmap")',
           call. = FALSE)
    } else stop('Requires R package "', package, '". Please install it.',
                call. = FALSE)
  }
}

check_chr_names <- function(gr) {
  check_package("GenomicRanges")
  if (!is(gr, "GRanges")) stop("Input must be a GRanges object")
 
  # Check only first few rows to speed up runtime 
  if (any(grep('chr[XVI]',
               as.character(GenomicRanges::seqnames(gr[1:10]))))) {
    return('roman numerals')
  } else if (any(grep('chr[0-9]',
                      as.character(GenomicRanges::seqnames(gr[1:50]))))) {
    return('arabic numerals')
  } else stop ('Did not recognize chromosome numbering system\n',
               'Please ensure chromosome numbers are in the usual format:\n',
               '"chrI" or "chr01".')
}

check_path <- function(path) {
  if (!file.exists(path)) {
    stop('Cannot seem to find file:\n', path,
         '\nPlease check that the provided path is correct.', call. = FALSE)
  }
}

make_local_copy <- function(path) {
  message('(copying file(s) to local folder "/hwglabr2_imports_temp"...)')
  # Check if directory already exists
  if (file.exists('hwglabr2_imports_temp')) {
    stop('A folder called "hwglabr2_imports_temp" already exists in the current',
         ' working directory.\n',
         'Please remove it and repeat function call.', call. = FALSE)
  }
  # Create temporary directory in cwd and make it the destination
  dir.create('hwglabr2_imports_temp')
  # Copy the files to the new temporary directory
  file.copy(path, 'hwglabr2_imports_temp', recursive = TRUE)
  # Update path to be the local directory
  path <- paste0('hwglabr2_imports_temp/', list.files('hwglabr2_imports_temp'))
  
  return(path)
}

check_gr_column <- function(gr, pattern='Y[A-P][LR]') {
  check_package("GenomicRanges")
  if (!is(gr, "GRanges")) stop("Input must be a GRanges object")
  
  metadata_columns <- names(GenomicRanges::mcols(gr))
  columns <- vector()
  for (column in metadata_columns) {
    if (any(grep(pattern,
                 as.character(as.data.frame(GenomicRanges::mcols(gr)[column]))))) {
      columns <- c(columns, column)
    }
  }
  return(columns)
}
