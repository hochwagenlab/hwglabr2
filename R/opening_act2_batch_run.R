#' Batch run of \code{\link{opening_act2}} function
#'
#' This function will take a CSV file containing the path to multiple data
#' folders and argument labels and run \code{\link{opening_act2}} function for
#' each data set in the file. It allows you to run a batch of data sets
#' automatically.
#' @param input_data_file String indicating path to file containing
#' instructions. The file must be in CSV format with the following columns
#' (column names in the file must be reproduced exactly as below):
#' \enumerate{
#'   \item \strong{path} Path to bedGraph file.
#'   \item \strong{genotype} String indicating the relevant strain
#'   mutations, used as the namesake argument to \code{\link{opening_act2}}.
#'   \item \strong{chip_target} String indicating the ChIP target protein, to be
#'   used as the namesake argument to \code{\link{opening_act2}}.
#'   \item \strong{sample_id} The sample ID, to be used as the namesake argument
#'   to \code{\link{opening_act2}}.
#' }
#' No default.
#' @param ref_genome Character object specifying the genome version, to be used
#' as the \code{genome} argument to \code{\link{opening_act2}}; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#' }
#' No default.
#' @param output_path Character object with a valid path to directory to save
#' output files at. Defaults to
#' \code{'/Volumes/LabShare/HTGenomics/Opening_act/'}.
#' @examples
#' \dontrun{
#' opening_act_batch_run(input_data_file="~/Desktop/inputs.csv",
#'                       ref_genome="SK1Yue")
#' 
#' opening_act_batch_run(input_data_file="~/Desktop/inputs.csv",
#'                       ref_genome="sacCer3", output_path='~/Desktop/')
#' }
#' @export

opening_act2_batch_run <- function(input_data_file, ref_genome,
                                   output_path=
                                     '/Volumes/LabShare/HTGenomics/Opening_act/'){
  t0 <- proc.time()[3]
  
  # IO checks
  check_path(input_data_file)
  
  # Read in file with path to every wiggle data set and all annotations
  args_file <- read.csv(input_data_file, stringsAsFactors = F)
  
  # Are all "paths" accessible?
  if (any(!file.exists(args_file[, 'path']))) {
    stop('Cannot seem to find the files for sample(s): ',
         paste(which(!file.exists(args_file[, 'path'])), collapse=", "), '\n',
         '       Please check the provided path(s) to files.', call. = FALSE)
  }
  
  # Are all "paths" unique?
  if (length(unique(args_file[, 'path'])) != nrow(args_file)) {
    stop('Not all provided paths to data are unique.\n',
         '       Please check the provided paths.',
         call. = FALSE)
  }
  
  # Have users check the arguments file for mistakes:
  # Ask user to make sure they provided valid arguments for opening_act
  title <- paste0('The information in the CSV file will be used to name the',
                  ' final output folders. For example, the "sampleID" label',
                  ' should identify the yeast strain, date, and read mapping',
                  ' conditions, as in:\n\n', '    AH7797-032817-SK1Yue-PM\n\n',
                  'Please make sure all your provided labels are correct.')
  choices <- c('Need to make changes, please quit.',
              'Looking good, continue analysis!')
  answer <- menu(choices, graphics = FALSE, title)
  
  if(answer == 0 | answer == 1){
    stop('You chose to quit. See you soon!', call. = FALSE)
  }
  
  # For each data set, read in wiggle data and run 'opening_act2'
  for(i in 1:nrow(args_file)){
    data <- hwglabr2::import_bedGraph(path=args_file[i, 'path'],
                                      local_copy=TRUE, keep_zeros=FALSE)
    hwglabr2::opening_act2(signal_data=data, genome=ref_genome,
                           genotype=args_file[i, 'genotype'],
                           chip_target=args_file[i, 'chip_target'],
                           sample_id=args_file[i, 'sample_id'],
                           output_path=output_path, user_input=FALSE)
    
    message('---')
    message('---')
    message('Ran "opening_act2" on all data found in ', args_file[i, 'path'])
  }
  
  message()
  message('----------------------------------------')
  message('----------------------------------------')
  message('Ran "opening_act" on all ', nrow(args_file), ' datasets.')
  message()
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
  message('----------------------------------------')
}