#------------------------------------------------------------------------------#
#                                                                              #
#                  Generation of data internal to the package                  #
#                                                                              #
#------------------------------------------------------------------------------#

# Get all data frames and then generate internal package data at the end of the
# script (/R/sysdata.rda)


#------------------------------------------------------------------------------#
#                               Centromeres                                    #

# SK1 genome assembly published in Yue et al. 2017
# Chromosome lengths calculated using:
#cat SK1.genome.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } \
#$0 !~ ">" {c+=length($0);} END { print c; }'
SK1Yuecen <- data.frame("Chromosome" = c("chrI","chrII","chrIII","chrIV","chrV",
                                         "chrVI","chrVII","chrVIII","chrIX",
                                         "chrX","chrXI","chrXII","chrXIII",
                                         "chrXIV","chrXV","chrXVI"),
                     "Start" = c(154628, 251815, 108708, 460752, 171136, 170910,
                                 501251, 102251, 348027, 447447, 451859, 151679,
                                 251000, 637019, 307189, 555578),
                     "End" = c(154745, 251931, 108824, 460871, 171253, 171027,
                               501370, 102368, 348143, 447565, 451975, 151797,
                               251118, 637136, 307307, 555694),
                     "LenChr" = c(228861, 829469, 340914, 1486921, 589812,
                                  299318, 1080440, 542723, 449612, 753937,
                                  690901, 1054145, 923535, 791982, 1053869,
                                  946846))

SK1Yuecen <- with(SK1Yuecen,
                  GenomicRanges::GRanges(Chromosome,
                                         IRanges::IRanges(Start + 1, End),
                                         seqlengths=setNames(LenChr,
                                                             Chromosome)))

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
#                         Intergenic region data                               #
#                         (Conv, div and tandem)                               #
# The intergenic region coordinate file generated using the scripts in:
# '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr2/'
# S288C and old SK1 regions were generated from the files in hwglabr

# 1. Import SK1Yue data
path <- '/Volumes/LabShare/GenomeSequences/hwglabr2/'
SK1Yue_intergenic <- read.table(paste0(path, 'SK1Yue_intergenic.txt'),
                                header = TRUE, stringsAsFactors = FALSE)
# 2. Import SK1 data
SK1_intergenic <- read.table(paste0(path, 'SK1_intergenic.txt'),
                             header = TRUE, stringsAsFactors = FALSE)
# 3. Import S288C data
sacCer3_intergenic <- read.table(paste0(path, 'sacCer3_intergenic.txt'),
                                 header = TRUE, stringsAsFactors = FALSE)

#------------------------------------------------------------------------------#
#                           Red1 summits in WT                                 #
#                         (Conv, div and tandem)                               #
# Files generated using the scripts in:
# '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr2/'
# 1. Import SK1Yue data
path <- '/Volumes/LabShare/GenomeSequences/hwglabr2/'
SK1Yue_Red1_summits_file <- paste0(path,
                                   'Red1-wildtype-71-34-199-29-Reps-SacCer3-',
                                   '2mis_B3W3_MACS2_over20_summits.bed')
SK1Yue_Red1_summits <- rtracklayer::import.bed(SK1Yue_Red1_summits_file)

# 2. Import S288C data
sacCer3_Red1_summits_file <- paste0(path,
                                    'Red1-wildtype-71-34-199-29-Reps-SK1Yue-',
                                    'PM_B3W3_MACS2_over20_summits.bed')
sacCer3_Red1_summits <- rtracklayer::import.bed(sacCer3_Red1_summits_file)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                           Add all data to package                            #
#                             (as internal data)                               #

# Determine the best compression for the data files
tools::checkRdaFiles('R/') # Suggests 'bzip2'

# Set package directory as working directory
# setwd('/path/to/hwglabr2/')
devtools::use_data(SK1Yuecen, sacCer3cen, SK1cen,
                   SK1Yue_intergenic, SK1_intergenic, sacCer3_intergenic,
                   SK1Yue_Red1_summits, sacCer3_Red1_summits,
                   internal = TRUE, overwrite = TRUE, compress = "bzip2")
