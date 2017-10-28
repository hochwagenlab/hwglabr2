#' FACS density plot
#'
#' Given a path to a directory containing \code{.fcs} files and a pattern to
#' match in the file names, generates density plots of FACS data across time
#' points.
#' @param dir Path to directory containing \code{.fcs} files. No default.
#' @param strain_code The lab's database numeric code identifier of the yeast
#' strain. No default.
#' @param gate Optional two-element vector containing gating limits. No default
#' (corresponds to no gating).
#' @param plot_color Optional R color to fill the plot. Defaults to
#' \code{black}.
#' @param plot_transparency Floating point number between 0 and 1 indicating 
#' transparency of the plot (corresponds to the plotting \code{alpha} argument:
#' 0 is full transparency and 1 is completely opaque). Defaults to \code{0.8}.
#' (corresponds to no gating).
#' @param file_format Character string indicating the plot file format. Must be
#' one of the options \code{jpeg} and \code{pdf}. Defaults to \code{jpeg}.
#' @param user_input Logical indicating whether to ask user before saving plot
#' to file. Defaults to \code{TRUE}.
#' @section Details:
#' The function is designed to search all \code{.fcs} files (named according to
#' the convention "Well strain_timepoint.fcs". e.g. "A01 119_0.fcs") for those
#' matching the provided pattern (yeast strain code) and build a density plot.
#' An ungated version of the plot is printed to the screen and the final plot is
#' saved to a file in the selected format in the same directory as the
#' \code{.fcs} files.
#' @return Density plot file in the same directory as the \code{.fcs} files.
#' @examples
#' \dontrun{
#' facs_density_plot2(119)
#' 
#' facs_density_plot2(dir='~/Desktop', strain_code=119,
#'                    gate=c(1500000, 7000000), file_format='pdf')
#'
#' facs_density_plot2(dir='~/Desktop', strain_code=119,
#'                    gate=c(1500000, 7000000), plot_color='orange',
#'                    plot_transparency=1)
#' }
#' @export

facs_density_plot2 <- function(dir, strain_code, gate, plot_color='black',
                               plot_transparency=0.8, file_format='jpeg',
                               user_input=TRUE){
  
  # IO checks
  if (is(dir, "character") & length(list.files(dir)) > 0) {
    check_path(dir)
  } else stop("'path' argument must be a path to a bedGraph file")
  
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop('Requires R package "flowCore". Please install it.\n',
         '          source("http://bioconductor.org/biocLite.R")\n',
         '          biocLite("flowCore")', call. = FALSE)
  }
  
  if (!requireNamespace("ggcyto", quietly = TRUE)) {
    stop('Requires R package "ggcyto". Please install it.\n',
         '          source("http://bioconductor.org/biocLite.R")\n',
         '          biocLite("ggcyto")', call. = FALSE)
  }
  
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop('Requires R package "ggridges". Please install it.\n',
         '          install.packages("ggridges")', call. = FALSE)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('Requires R package "ggplot2". Please install it.\n',
         '          install.packages("ggplot2")', call. = FALSE)
  }
  
  if(!file_format %in% c('jpeg', 'pdf')){
    stop('"file_format" must be one of "jpg" and "pdf".', call. = FALSE)
  }
  
  library("flowCore")
  library("ggcyto")
  library("ggridges")
  #library("ggplot2")
  
  # Import .fcs files
  message('Loading fcs files...')
  files <- list.files(dir)
  target_files <- files[grep(strain_code, files)]
  target_files <- target_files[grep("fcs", target_files)]
  
  # Load all files into list
  samples <- vector("list")
  for (file in target_files) {
    # Determine the time point of the file. This assumes that the format of 
    # the file name is: date_yeastLine_timePoint.fcs
    time_point <- strsplit(strsplit(file, "\\.") [[1]][1], "_")[[1]][2]
    
    # Read the FCS file
    samples[[time_point]] <- flowCore::read.FCS(paste0(dir, file))
  }
  
  # Convert to "flowSet"
  fs <- as(samples, "flowSet")
  
  message('Plotting data...')
  # Plot
  p <- ggplot(fs, aes(x=`FL1-A`, y=name)) +
    geom_density_ridges(fill=plot_color, scale=2, color='white',
                        size=0.25, alpha=plot_transparency) +
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_fill_gradient2() +
    labs(title = 'WT', x = "DNA content", y = "Time in meiosis\n(hrs)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.4, size = 30, face = 'bold'),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = 'black', size = 26),
          axis.title = element_text(size = 26))
  
  # Show ungated plot on screen
  print(p)
  
  # Ask user to save file
  if (missing(gate)) {
    question <- 'Save plot to file?'
  } else question <- 'Gate plot and save to file?'
  
  if(user_input) {
    title <- question
    choices = c('No, let me change it first.', 'Yes, save plot!')
    answer <- menu(choices, graphics = FALSE, title)
    
    if(answer == 0 | answer == 1){
      stop('You chose to stop the function.', call. = FALSE)
    }
  }
  
  file_name <- paste0(dir, strain_code)
  
  if (!missing(gate)) print(p + xlim(gate))
  message('Saving plot:')
  
  if (file_format == 'jpeg') {
    message('   ', paste0(file_name, '.jpg'))
    dev.copy(jpeg, paste0(file_name, '.jpg'))
  } else {
    message('   ', paste0(file_name, '.pdf'))
    dev.copy(pdf, paste0(file_name, '.pdf'))
  }
  
  dev.off()
  
  message('---')
  message('Done!')
}