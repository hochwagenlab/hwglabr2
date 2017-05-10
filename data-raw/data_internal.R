#------------------------------------------------------------------------------#
#                                                                              #
#                  Generation of data internal to the package                  #
#                                                                              #
#------------------------------------------------------------------------------#

# Get all data frames and then generate internal package data at the end of the
# script (/R/sysdata.rda)


#------------------------------------------------------------------------------#
#                               Centromeres                                    #

# SK1 info based on Keeney lab genome sequence and annotation
SK1cen <- data.frame("Chromosome" = c("chr01","chr02","chr03","chr04","chr05",
                                     "chr06","chr07","chr08","chr09","chr10",
                                     "chr11","chr12","chr13","chr14","chr15",
                                     "chr16"),
                     "Start" = c(137832, 226711, 128699, 463204, 157003, 162815,
                                 505440, 95031, 346215, 415648, 452723, 137738,
                                 249103, 616840, 307236, 553355),
                     "End" = c(137948, 226826, 128779, 463321, 157119, 162931,
                               505558, 95147, 346330, 415764, 452838, 137855,
                               249221, 616956, 307353, 553467),
                     "Mid" = c(137890, 226768, 128739, 463262, 157061, 162873,
                               505499, 95089, 346272, 415706, 452780, 137796,
                               249162, 616898, 307294, 553411),
                     "LenChr" = c(203893, 794508, 342718, 1490682, 602514,
                                  284456, 1067526, 544538, 435585, 719294,
                                  687260, 1008248, 908607, 812465, 1054033,
                                  921188))

SK1cen <- with(SK1cen, GenomicRanges::GRanges(Chromosome,
                                              IRanges::IRanges(Start + 1, End),
                                              seqlengths=setNames(LenChr,
                                                                  Chromosome)))

# S288c
path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/'
file <- 'saccharomyces_cerevisiae_R64-1-1_20110208.gff'
sacCer3gff <- read.table(paste0(path, file),
                         fill = TRUE, stringsAsFactors = FALSE)
sacCer3gff <- sacCer3gff[1:16406, ]
sacCer3cen <- sacCer3gff[sacCer3gff[, 3] == 'centromere', c(1, 4:5)]
names(sacCer3cen) <- c('Chromosome', 'Start', 'End')
sacCer3cen$Mid <- floor(sacCer3cen$Start + (sacCer3cen$End - sacCer3cen$Start) / 2)
sacCer3cen$LenChr <- sacCer3gff[sacCer3gff[, 3] == 'chromosome', 5][1:16]

sacCer3cen <- with(sacCer3cen,
                   GenomicRanges::GRanges(Chromosome,
                                          IRanges::IRanges(Start+1, End),
                                          seqlengths = setNames(LenChr,
                                                                Chromosome)))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                           Add all data to package                            #
#                             (as internal data)                               #

# Determine the best compression for the data files
tools::checkRdaFiles('R/') # Suggests 'bzip2'

# Set package directory as working directory
# setwd('/path/to/hwglabr2/')
devtools::use_data(sacCer3cen, SK1cen,
                   internal = TRUE, overwrite = TRUE, compress = "bzip2")
