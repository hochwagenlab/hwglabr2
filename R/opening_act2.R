#' Standard analysis of ChIP-seq experiments
#'
#' Runs the lab's standard analysis of ChIP-seq experiment data and produces
#' several .pdf files of analysis plots to be kept in a new folder at
#' \code{.../LabShare/HTGenomics/Opening_act/}.
#' 
#' @param signal_data Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1Yue"}
#'   \item \code{"sacCer3"}
#' }
#' No default.
#' @param genotype Character object indicating the relevant strain mutations.
#' Just use \code{genotype=WT}, for example, if there are no relevant mutations.
#' No default.
#' @param chip_target Character object indicating the ChIP target protein.
#' No default.
#' @param sample_id Character object indicating the sample ID, including the ID
#' used in the analysis pipeline (with a date) and the read mapping conditions
#' (see examples below). If \code{user_input=TRUE} the function asks the user to
#' check that the provided \code{sample_id} matches the required format before
#' proceeding with the analysis. No default.
#' @param output_path Character object with a valid path to directory to save
#' output files at. Defaults to
#' \code{'/Volumes/LabShare/HTGenomics/Opening_act/'}.
#' Defaults to \code{output_path='/Volumes/LabShare/HTGenomics/Opening_act/'}
#' @param user_input Logical indicating whether to ask user to check the format
#' of the \code{sample_id} argument. Defaults to \code{TRUE}.
#' @param run_chr_size_bias Logical indicating whether to run the chromosome
#' size bias analysis. All included analysis can be turned off in this way.
#' Defaults to \code{TRUE}.
#' @param run_centromeres Logical indicating whether to run the signal at 
#' centromere analysis. Defaults to \code{TRUE}.
#' @param run_rDNA Logical indicating whether to run the signal flanking rDNA
#' analysis. Defaults to \code{TRUE}.
#' @param run_telomeres Logical indicating whether to run the signal at
#' sub-telomeric region analysis. Defaults to \code{TRUE}.
#' @param run_dsb_hotspots Logical indicating whether to run the signal at DSB
#' hotspot analysis. Defaults to \code{TRUE}.
#' @param run_axis Logical indicating whether to run the signal at axis binding
#' site analysis. Defaults to \code{TRUE}.
#' @param run_meta_orf Logical indicating whether to run the meta ORF analysis.
#' Defaults to \code{TRUE}.
#' 
#' @return A new folder at the specified location (default is
#' \code{.../LabShare/HTGenomics/Opening_act/}) containing output plots (as .pdf
#' files) of the included of the following analyses:
#' \enumerate{
#'   \item Chromosome size bias
#'   \item Signal at centromeres
#'   \item Signal flanking rDNA
#'   \item Signal at sub-telomeric regions
#'   \item Signal at DSB hotspots
#'   \item Signal at axis binding sites
#'   \item Signal at meta ORF
#' }
#' @examples
#' \dontrun{
#' opening_act2(signal_data=WT, genome="sacCer3", genotype="WT",
#'              chip_target="Red1", sample_id="AH119C-040114-sacCer3-2mis")
#'
#' opening_act2(set1_wiggle, genome="sacCer3", "set1", "Red1",
#'              "AH8584b-16032016-sacCer3-2mis")
#' 
#' opening_act2(rec8, genome="SK1Yue", "rec8", "Red1",
#'              "AH8115b-24042015-SK1Yue-PM", user_input=FALSE,
#'              run_meta_orf=FALSE)
#'             
#' opening_act2(signal_data=WT, genome="sacCer3", genotype="WT",
#'              chip_target="Red1", sample_id="AH119C-040114-sacCer3-2mis",
#'              run_meta_orf=FALSE)
#'
#' opening_act2(signal_data=WT, genome="sacCer3", genotype="dot1",
#'              chip_target="Red1", sample_id="AH8104-010116-sacCer3-2mis",
#'              output_path='~/Desktop',
#'              run_chr_size_bias=FALSE, run_centromeres=FALSE, run_rDNA=FALSE,
#'              run_telomeres=FALSE, run_dsb_hotspots=TRUE, run_axis=TRUE,
#'              run_meta_orf=FALSE)
#' }
#' @export

opening_act2 <- function(signal_data, genome, genotype, chip_target, sample_id,
                         output_path='/Volumes/LabShare/HTGenomics/Opening_act/',
                         user_input=TRUE,
                         run_chr_size_bias=TRUE, run_centromeres=TRUE,
                         run_rDNA=TRUE, run_telomeres=TRUE,
                         run_dsb_hotspots=TRUE, run_axis=TRUE,
                         run_meta_orf=TRUE) {
  t0 <- proc.time()[3]
  
  # IO checks
  check_package("GenomicRanges")
  check_package("EnrichedHeatmap")
  
  check_path(output_path)
  
  if (!is(signal_data, "GRanges")) {
      stop('"signal_data" must be a GRanges object.', call. = FALSE)
  }
  
  if (!genome %in% c('SK1Yue', 'sacCer3')) {
    stop('"genome" must be either "SK1Yue" or "sacCer3".', call. = FALSE)
  }
  
  # Ask user to double-check "sample_id" and "genome" argument input
  if(user_input){
    # Ask user to make sure they provided a valid ID for the data set
    title <- paste0('The "sampleID" argument will be used to name the final',
                    ' output folder.\n',
                    'It should identify the yeast strain, date, and read ',
                    'mapping conditions, as in:\n',
                    '"AH119C-040114-sacCer3-2mis"\n',
                    'You provided the following "sample_id" and "genome":\n',
                    '     ', sample_id, '\n     ', genome,
                    '\n\nIs this correct?')
    choices = c('No, let me change that.', 'Yes, continue analysis!')
    answer <- menu(choices, graphics = FALSE, title)
    
    if(answer == 0 | answer == 1){
      stop('You chose to stop the function.', call. = FALSE)
    }
  }
  
  # Create output directory (if it doesn't already exist)
  output_dir <- paste0(genotype, '_anti-', chip_target, '_', sample_id)
  if (file.exists(paste0(output_path, output_dir))) {
      stop('A folder named "', output_dir, '" already exists at\n', output_path,
           call. = FALSE)
  }
  message('Creating output directory "', output_dir, '"...')
  dir.create(file.path(paste0(output_path, output_dir)))
  
  #----------------------------------------------------------------------------#
  #                                Run analysis                                #
  #----------------------------------------------------------------------------#
  message('Running analyses:')
  
  #----------------------------------------------------------------------------#
  # Chr size bias
  
  if (run_chr_size_bias) {
    message('   Chromosome size bias')
    suppressMessages(output <- hwglabr2::average_chr_signal(signal_data,
                                                            remove_cen=F,
                                                            mean_norm=F)[[1]])
    # Get chr lengths and add to average signal
    chr_lengths <- GenomeInfoDb::seqlengths(
      hwglabr2::get_chr_coordinates(genome=genome)
      )
    chr_lengths <- data.frame(chr=names(chr_lengths), length=chr_lengths)
    output <- merge(output, chr_lengths)
    
    # Get genome-wide mean
    genome_wide_mean <- suppressMessages(
      hwglabr2::average_chr_signal(signal_data)[[2]]
    )
    
    # Normalize signal
    output$avrg_signal <- output$avrg_signal / genome_wide_mean
    
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_chrSizeBias.pdf')
    pdf(file=file_name, width=4, height=4)
    par(mar = c(5, 5, 4, 2))
    plot(output$length / 1000, output$avrg_signal,
         xlim=c(0, max(output$length) / 1000),
         xlab = 'Chromosome size (kb)', ylab = 'Signal (genome-wide mean = 1)',
         main = paste0('Mean signal / chromosome\nrelative to chromosome size'),
         col = 'grey50', pch = 19)
    abline(h = 1, lty=3, lwd=1.5)
    dev.off()
    
    message('      Saved ', paste0(output_dir, '_chrSizeBias.pdf'))
  } else {
    message('      (Skip chromosome size bias)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal at centromeres
  if (run_centromeres) {
    message('   Signal at centromeres')
    
    # Get centromeres
    centromeres <- hwglabr2::get_chr_coordinates(genome=genome, as_df=FALSE)
    
    # Replace start and end centromere positions by midpoint
    midpoint <- floor(GenomicRanges::width(centromeres) / 2)
    GenomicRanges::start(centromeres) <- GenomicRanges::start(centromeres) +
      midpoint
    GenomicRanges::end(centromeres) <- GenomicRanges::start(centromeres)
    
    signal_at_cens <- EnrichedHeatmap::normalizeToMatrix(signal_data,
                                                         centromeres,
                                                         extend=40000, w=50,
                                                         mean_mode="weighted",
                                                         value_column="score")
    
    signal_at_cens_avrg <- colMeans(signal_at_cens, na.rm = T)
    
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_signalAtCen.pdf')
    pdf(file = paste0(file_name), width = 6, height = 3)
    ylim <- range(signal_at_cens_avrg)
    if( ylim[2] < 2) ylim[2] <- 2
    plot(seq(-39999, 40000, 50) / 1000, signal_at_cens_avrg, type="l",
         ylim=ylim, xlab="Distance to centromere (kb)", ylab="Signal", 
         lwd=1, cex.axis=1, las=1, col="darkorange", cex.lab=1,
         main='Average signal around centromeres', cex.main=1)
    
    dev.off()
    message('      Saved ', paste0(output_dir, '_signalAtCen.pdf')) 
  } else {
    message('      (Skip signal at centromeres)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal at rDNA
  if (run_rDNA) {
    message('   Signal flanking rDNA')
    
    rDNA <- suppressMessages(hwglabr2::signal_flanking_rDNA(signal_data,
                                                            flank_length=40000,
                                                            genome=genome))
    
    # Set chromosome length to the last position in the collected data
    # (to avoid creating genome tiles for positions downstream of the rDNA)
    last_position <- tail(GenomicRanges::end(rDNA), 1)
    seq_length <- c('chrXII'=last_position)
    
    # Compute 50-bp tiling windows (will compress the data a little bit)
    bins <- GenomicRanges::tileGenome(seq_length, tilewidth=40,
                                      cut.last.tile.in.chrom=TRUE)
    
    # Keep only required region
    first_position <- head(GenomicRanges::start(rDNA), 1)
    greater_than_start <- GenomicRanges::start(bins) >= first_position
    smaller_than_end <- GenomicRanges::end(bins) <= last_position
    bins <- bins[greater_than_start & smaller_than_end]
    # Keep signal data for chr XII only (object's `seqlevels` must match bins')
    rDNA <- GenomeInfoDb::keepSeqlevels(rDNA, "chrXII")
    # Get signal as "RleList"; the signal is in the "score" metadata column
    score <- GenomicRanges::coverage(rDNA, weight="score")
    # Compute average signal per tile
    bins <- GenomicRanges::binnedAverage(bins, score, "binned_score")
    # Make data frame (convert positions to Kb; signal is the binned score)
    rDNA <- data.frame(position=GenomicRanges::start(bins),
                       signal=bins$binned_score)
    
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_signalAtrDNA.pdf')
    pdf(file = paste0(file_name), width = 6, height = 3)
    plot(rDNA$position / 1000, rDNA$signal, type="l",
         xlab="Position on chr XII (kb)", ylab="Signal", 
         lwd=1, cex.axis=1, las=1, col='black', cex.lab=1, cex.main=1)
    
    # A stretch present in S288C downstream of rDNA is absent in SK1
    # Add label for that (from end of rDNA until about bp 490'500)
    if (genome == 'sacCer3') {
      start <- 468931
      end <- 490500
      axis(1, at = c(start / 1000, end / 1000),
           labels = c('', ''),
           col = 'blue', lwd = 2)
    }
    
    # Add labels for rDNA
    if (genome == 'sacCer3') {
      start <- 451575
      end <- 468931
      title(paste0("Signal around rDNA: ",
                   '\n(rDNA position marked in red;',
                   '\nregion absent form SK1 genome marked in blue)'))
    } else {
      start <- 447012
      end <- 461699
      title(paste0("Signal around rDNA: ", '\n(rDNA position marked in red)'))
    }
    
    axis(1, at = c(start / 1000, end / 1000),
         labels = c('', ''),
         col = 'red', lwd = 2)
    
    dev.off()
    message('      Saved ', paste0(output_dir, '_signalAtrDNA.pdf'))
  } else {
    message('      (Skip signal flanking rDNA)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal from telomeres
  if (run_telomeres) {
    message('   Signal at sub-telomeric regions')
    
    signal_at_tels <- suppressMessages(
      hwglabr2::signal_from_telomeres2(signal_data, 120000, 500, genome)
      )
    
    signal_at_tels <- colMeans(signal_at_tels[, 4:ncol(signal_at_tels)],
                               na.rm = T)
    
    # Get genome-wide mean
    genome_wide_mean <- suppressMessages(
      hwglabr2::average_chr_signal(signal_data)[[2]]
      )
    
    # Normalize signal
    signal_at_tels <- signal_at_tels / genome_wide_mean
    
    # Plot results
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_signalAtTelomeres.pdf')
    pdf(file = paste0(file_name), width = 6, height = 4)
    
    plot(x=seq(1, 120000, 500), signal_at_tels, type='l', lwd=2, col='plum4',
         xlab='Distance from chr end (Kb)', ylab='Signal (genome-wide mean = 1)',
         main='Signal at sub-telomeric regions',
         cex.main=1)
    abline(h = 1, lty=3, lwd=1.5)
    dev.off()
    message('      Saved ', paste0(output_dir, '_signalAtTelomeres.pdf'))
  } else {
    message('      (Skip signal at sub-telomeric regions)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal around DSBs by DSB hotspot hotness
  if (run_dsb_hotspots) {
    message('   Signal at DSB hotspots')
    
    Spo11_hs <- hwglabr2::get_dsb_hotspots(genome)
    midpoint <- floor(GenomicRanges::width(Spo11_hs) / 2)
    GenomicRanges::start(Spo11_hs) <- GenomicRanges::start(Spo11_hs) + midpoint
    GenomicRanges::end(Spo11_hs) <- GenomicRanges::start(Spo11_hs)
    
    # Order by signal to make 8 groups based on hotspot hotness
    Spo11_hs <- Spo11_hs[order(Spo11_hs$score), ]
    
    signal_at_hs <- EnrichedHeatmap::normalizeToMatrix(signal_data, Spo11_hs,
                                                       extend=1000, w=10,
                                                       mean_mode="weighted",
                                                       value_column="score")
    
    # Define quantiles (matrix row order respects order by score from above)
    hs_quants <- list()
    number_of_hs <- nrow(signal_at_hs)
    quant_length <- round(number_of_hs / 8)
    for (i in 1:8) {
      start <- 1 + ((i - 1) * quant_length)
      end <- min(number_of_hs, (start + quant_length) - 1)
      hs_quants[[i]] <- signal_at_hs[start:end, ]
    }
    
    # Average data
    hs_quants <- lapply(hs_quants, function(x) colMeans(x, na.rm = T))
    
    min_data <- sapply(hs_quants, function(x) min(x))
    max_data <- sapply(hs_quants, function(x) max(x))
    
    colors <- c('lightblue1', 'cadetblue1', 'deepskyblue', 'deepskyblue3',
                'royalblue', 'blue', 'blue4', 'black')
    
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_signalAtDSBhotspots.pdf')
    pdf(file = paste0(file_name), width = 6, height = 5)
    
    plot(0, type='l', lwd=3, xlim = c(-1000, 1000), ylab='Signal',
         ylim=c(min(min_data), max(max_data)),
         xlab='Distance from DSB hotspot midpoints (bp)',
         main = 'Signal around DSB hotspots')
    
    for(i in 1:length(hs_quants)){
      lines(x=seq(-999, 1000, 10), y=hs_quants[[i]], lwd=3, col=colors[i])
    }
    
    legend('topright', lty=c(1,1), lwd=3, title='Hotspot strength',
           legend=c('weakest', '', '', '', '', '', '', 'hottest'),
           bg = 'white', col=colors, cex=0.6)
    dev.off()
    
    message('      Saved ', paste0(output_dir, '_signalAtDSBhotspots.pdf'))
  } else {
    message('      (Skip signal at DSB hotspots)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal at axis binding sites
  if (run_axis) {
    message('   Signal at axis binding sites')
    
    Red1_summits <- hwglabr2::get_Red1_summits(genome)
    # Order by signal to make groups based on binding level
    Red1_summits <- Red1_summits[order(Red1_summits$score), ]
    
    signal_at_axis <- EnrichedHeatmap::normalizeToMatrix(signal_data,
                                                         Red1_summits,
                                                         extend=1000, w=10,
                                                         mean_mode="weighted",
                                                         value_column="score")
    
    # Define quantiles (matrix row order respects order by score from above)
    axis_quants <- list()
    number_of_peaks <- nrow(signal_at_axis)
    quant_length <- round(number_of_peaks / 4)
    for (i in 1:4) {
      start <- 1 + ((i - 1) * quant_length)
      end <- min(number_of_peaks, (start + quant_length) - 1)
      axis_quants[[i]] <- signal_at_axis[start:end, ]
    }
    
    # Average data
    axis_quants <- lapply(axis_quants, function(x) colMeans(x, na.rm = T))
    
    min_data <- sapply(axis_quants, function(x) min(x))
    max_data <- sapply(axis_quants, function(x) max(x))
    
    colors <- c("darkolivegreen3", "green3", "darkgreen", "black")
    
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_signalAtAxisSites.pdf')
    pdf(file = paste0(file_name), width = 6, height = 5)
    
    plot(0, type='l', lwd=3, xlim = c(-1000, 1000), ylab='Signal',
         ylim=c(min(min_data), max(max_data)),
         xlab='Distance from Axis Midpoints (bp)',
         main = 'Signal around axis sites\n(Red1 peaks in WT)')
    
    for(i in 1:length(axis_quants)){
      lines(x=seq(-999, 1000, 10), y=axis_quants[[i]], lwd=3, col=colors[i])
    }
    
    legend('topright', lty=c(1, 1), lwd=3, title='Axis binding',
           legend=c('least', '', '', 'most'), bg='white', col=colors, cex=0.6)
    dev.off()
    
    message('      Saved ', paste0(output_dir, '_signalAtAxisSites.pdf'))
    
  } else {
    message('      (Skip signal at axis binding sites)')
  }
  
  #----------------------------------------------------------------------------#
  # Meta ORF
  if (run_meta_orf) {
    message('   Signal at meta ORF')
    
    gff <- hwglabr2::get_gff(genome)
    
    # Drop unnecessary annotations (only actually applies with SK1Yue)
    gff <- gff[gff$type %in% c('gene', 'pseudogene')]
    
    signal_at_ORFs <- suppressMessages(
      hwglabr2::signal_at_orf2(signal_data, gff, write_to_file=FALSE)
      )
    
    meta_orf <- colMeans(signal_at_ORFs, na.rm=TRUE)
    
    # plot results
    file_name <- paste0(output_path, output_dir, '/', output_dir,
                        '_signalAtORF.pdf')
    pdf(file = paste0(file_name), width = 6, height = 4)
    plot(meta_orf, type='l', xaxt='n', yaxt='n',
         xlim=c(0, 1000), lwd = 3, col = 'orange',
         xlab="Scaled ORF",
         ylab='Signal', main=paste0('Signal at meta ORF'), bty='n')
    axis(1, at=c(0, 250, 750, 1000), labels=c('', 'start', 'stop', ''), las=1)
    axis(2, las=2)
    abline(v=c(250, 750), lty=2)
    dev.off()
    message('      Saved ', paste0(output_dir, '_signalAtmetaORF.pdf')) 
  } else {
    message('      (Skip signal at meta ORF)')
  }
  #----------------------------------------------------------------------------#
  
  message()
  message('------------------')
  message('All plots saved to ', paste0(output_path, output_dir))
  message('------------------')
  message()
  message('---')
  message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
}
