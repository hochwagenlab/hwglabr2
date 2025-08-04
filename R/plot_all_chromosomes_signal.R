#' Plot genome-wide signal across chromosomes
#'
#' Bins signal genome-wide and plots individual chromosome-level signal tracks 
#' in a vertical stack. This is useful for visualizing genome-wide patterns from 
#' ChIP-seq or related assays across all chromosomes in a fixed layout.
#'
#' @param sample_bdg A \code{GRanges} object containing genome-wide signal 
#' (e.g., from ChIP-seq). Must include a \code{score} metadata column and 
#' chromosome names prefixed with "chr".
#'
#' @param tile Integer. Genomic bin size (in bp) to average signal over. 
#' Defaults to \code{100}.
#'
#' @param ylim_vals Numeric vector of length 2. Fixed y-axis range for all 
#' panels. Defaults to \code{c(-0.65, 10)}.
#'
#' @param xlim_max Integer. Maximum x-axis (genomic coordinate) for all 
#' chromosomes. Used to standardize panel widths. Defaults to \code{1531933}.
#'
#' @param color Character string indicating the color of the centromere marker 
#' in each panel. Accepts hex codes or color names. Defaults to \code{"black"}.
#'
#' @return A \code{patchwork} object containing \code{ggplot2} plots for each 
#' chromosome, stacked vertically.
#'
#' @examples
#' \dontrun{
#' plot_all_chromosomes_signal(sample_bdg, tile = 100, color = "darkgreen")
#' }
#'
#' @export
plot_all_chromosomes_signal <- function(sample_bdg,
                                        tile = 100,
                                        ylim_vals = c(-0.65, 10),
                                        xlim_max = 1531933,
                                        color = "black") {
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ggplot2)
  library(patchwork)
  library(hwglabr2)

  # Set plotting theme
  ggplot2_blanktheme <- theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )
  theme_set(ggplot2_blanktheme)

  # Chromosome order (SK1Yue)
  chrs <- list('I', 'VI', 'III', 'IX', 'VIII', 'V', 'XI', 'X',
               'XIV', 'II', 'XIII', 'XVI', 'XII', 'VII', 'XV', 'IV')

  genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

  # Bin and plot helper
  plot_signal <- function(sample, chrnum, tile) {
    sample_ex <- sort(sortSeqlevels(sample))
    seqlengths(sample_ex) <- seqlengths(genome_info)
    bins_ex <- tileGenome(seqlengths(sample_ex),
                          tilewidth = tile,
                          cut.last.tile.in.chrom = TRUE)
    score_ex <- coverage(sample_ex, weight = "score")
    bins_ex <- binnedAverage(bins_ex, score_ex, "binned_score")
    bins_ex <- keepSeqlevels(bins_ex, paste0("chr", chrnum), pruning.mode = "coarse")
    positions_ex <- start(bins_ex) + floor(width(bins_ex) / 2)
    df_ex <- data.frame(
      seqnames = paste0("chr", chrnum),
      position = positions_ex,
      signal = bins_ex$binned_score
    )
    return(df_ex)
  }

  plist <- list()
  for (chr in chrs) {
    df <- plot_signal(sample_bdg, chr, tile)
    cen <- genome_info[seqnames(genome_info) == paste0("chr", chr)]
    cen_mid <- round((start(cen) + end(cen)) / 2)
    cen_df <- data.frame(cen_midpt = cen_mid, y = 0)

    g <- ggplot(df, aes(position, signal)) +
      geom_line() +
      ylab(chr) +
      xlim(0, xlim_max) +
      ylim(ylim_vals) +
      geom_point(
        data = cen_df,
        aes(x = cen_midpt, y = -0.15),
        size = 1.5,
        colour = color,
        shape = 20
      )

    plist[[chr]] <- g
  }

  return(wrap_plots(plist, nrow = length(plist)))
}
