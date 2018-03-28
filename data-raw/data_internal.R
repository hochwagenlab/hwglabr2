#------------------------------------------------------------------------------#
#                                                                              #
#                  Generation of data internal to the package                  #
#                                                                              #
#------------------------------------------------------------------------------#

# Get all data frames and then generate internal package data at the end of the
# script (/R/sysdata.rda)

#------------------------------------------------------------------------------#
#                             Helper functions                                 #
add_genome_name_to_GR <- function(gr, name='SK1Yue') {
  number_seqs <- length(levels(gr@seqnames))
  gr@seqinfo@genome <- rep(name, number_seqs)
  
  gr
}


#------------------------------------------------------------------------------#
#                                GFF files                                     #
path <- '/Volumes/LabShare/GenomeSequences'
# 1. Import SK1Yue data
SK1Yue_gff <- 'SK1_Yue_et_al_2017/Yue.SK1.genome.nuclear.mito.2micr.gff'
SK1Yue_gff <- rtracklayer::import.gff3(file.path(path, SK1Yue_gff))
SK1Yue_gff <- add_genome_name_to_GR(SK1Yue_gff, name='SK1Yue')
# 2. Import S288CYue data
# $ cd Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/
# $ mkdir S288C_Yue_et_al_2017/
# $ cd S288C_Yue_et_al_2017/
# $ wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/S288c.all_feature.gff.gz
# $ gunzip *
S288CYue_gff <- rtracklayer::import.gff(
  file.path('/Volumes/LabShare-2/Luis/LabWork/GenomeSequences',
            'S288C_Yue_et_al_2017/S288c.all_feature.gff'))
# 3. Import S288C data
sacCer3_gff <- 's288C_annotation_R64_modified.gff'
sacCer3_gff <- rtracklayer::import.gff(file.path(path, sacCer3_gff))
sacCer3_gff <- add_genome_name_to_GR(sacCer3_gff, name='sacCer3')

sacCer2_gff <- file.path(
  'S288C_reference_genome_R61-1-1_20080605',
  'saccharomyces_cerevisiae_R61-1-1_20080607_modified.gff')
sacCer2_gff <- rtracklayer::import.gff(file.path(path, sacCer2_gff))
sacCer2_gff <- add_genome_name_to_GR(sacCer2_gff, name='sacCer2')
# 4. Import SK1 data
SK1_gff <- 'SK1_MvO_V1___GENOME/SK1_annotation/SK1_annotation_modified_v2.gff'
SK1_gff <- rtracklayer::import.gff(file.path(path, SK1_gff))
SK1_gff <- add_genome_name_to_GR(SK1_gff, name='SK1')


#------------------------------------------------------------------------------#
#                               Centromeres                                    #

# SK1 genome assembly published in Yue et al. 2017
# Chromosome lengths calculated using:
#cat SK1.genome.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } \
#$0 !~ ">" {c+=length($0);} END { print c; }'
SK1Yuecen <- data.frame(
  "Chromosome" = paste0('chr', as.roman(1:16)),
  "Start" = c(154628, 251815, 108708, 460752, 171136, 170910, 501251, 102251,
              348027, 447447, 451859, 151679, 251000, 637019, 307189, 555578),
  "End" = c(154745, 251931, 108824, 460871, 171253, 171027, 501370, 102368,
            348143, 447565, 451975, 151797, 251118, 637136, 307307, 555694),
  "LenChr" = c(228861, 829469, 340914, 1486921, 589812, 299318, 1080440, 542723,
               449612, 753937, 690901, 1054145, 923535, 791982, 1053869,
               946846))

SK1Yuecen <- with(SK1Yuecen,
                  GenomicRanges::GRanges(Chromosome,
                                         IRanges::IRanges(Start, End),
                                         seqlengths=setNames(LenChr,
                                                             Chromosome)))

SK1Yuecen <- add_genome_name_to_GR(SK1Yuecen, name='SK1Yue')

# S288c genome assembly published in Yue et al. 2017
# Chromosome lengths calculated using:
#cat S288c_SK1_Yue.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } \
#$0 !~ ">" {c+=length($0);} END { print c; }'
start <- S288CYue_gff[S288CYue_gff$type == 'centromere']@ranges@start
end <- start + S288CYue_gff[S288CYue_gff$type == 'centromere']@ranges@width - 1

S288CYuecen <- data.frame(
  "Chromosome" = paste0('chr', as.roman(1:16)),
  "Start" = start, "End" = end,
  "LenChr" = c(219929, 813597, 341580, 1566853, 583092, 271539, 1091538, 581049,
               440036, 751611, 666862, 1075542, 930506, 777615, 1091343,
               954457))

S288CYuecen <- with(S288CYuecen,
                    GenomicRanges::GRanges(Chromosome,
                                           IRanges::IRanges(Start, End),
                                           seqlengths=setNames(LenChr,
                                                               Chromosome)))

S288CYuecen <- add_genome_name_to_GR(S288CYuecen, name='S288cYue')

# SK1 info based on Keeney lab genome sequence and annotation
SK1cen <- data.frame(
  "Chromosome" = paste0('chr', stringr::str_pad(1:16, 2, pad = "0")),
  "Start" = c(137832, 226711, 128699, 463204, 157003, 162815, 505440, 95031,
              346215, 415648, 452723, 137738, 249103, 616840, 307236, 553355),
  "End" = c(137948, 226826, 128779, 463321, 157119, 162931, 505558, 95147,
            346330, 415764, 452838, 137855, 249221, 616956, 307353, 553467),
  "Mid" = c(137890, 226768, 128739, 463262, 157061, 162873, 505499, 95089,
            346272, 415706, 452780, 137796, 249162, 616898, 307294, 553411),
  "LenChr" = c(203893, 794508, 342718, 1490682, 602514, 284456, 1067526, 544538,
               435585, 719294, 687260, 1008248, 908607, 812465, 1054033,
               921188))

SK1cen <- with(SK1cen, GenomicRanges::GRanges(Chromosome,
                                              IRanges::IRanges(Start, End),
                                              seqlengths=setNames(LenChr,
                                                                  Chromosome)))
SK1cen <- add_genome_name_to_GR(SK1cen, name='SK1')

# S288c
path <- file.path('/Volumes/LabShare/GenomeSequences',
                  'S288C_reference_genome_R64-1-1_20110203',
                  'saccharomyces_cerevisiae_R64-1-1_20110208.gff')
sacCer3gff <- rtracklayer::import.gff(path)
sacCer3cen <- sacCer3gff[sacCer3gff$type == 'centromere', ]
sacCer3cen <- GenomeInfoDb::dropSeqlevels(sacCer3cen, c('chrMito', '2-micron'))
GenomicRanges::mcols(sacCer3cen) <- NULL
chr_len <- sacCer3gff[sacCer3gff$type == 'chromosome'][1:16]
chr_len <- setNames(chr_len@ranges@width, chr_len$Name)
GenomeInfoDb::seqlengths(sacCer3cen) <- chr_len

sacCer3cen <- add_genome_name_to_GR(sacCer3cen, name='sacCer3')

### SK1-S288C hybrid genome using assemblies published in Yue et al. 2017
# Using the data from the SK1 and the S288C Yue et al. assemblies from above
# Useful for analysis of hybrid strains and SK1 strain ChIP-seq spiked in with
# SK288c cells

start <- S288CYue_gff[S288CYue_gff$type == 'centromere']@ranges@start
end <- start + S288CYue_gff[S288CYue_gff$type == 'centromere']@ranges@width - 1

SK1_S288CYuecen <- data.frame(
  "Chromosome" = c(paste0('chr', as.roman(1:16), '_SK1'),
                   paste0('chr', as.roman(1:16), '_S288C')),
  "Start" = c(c(154628, 251815, 108708, 460752, 171136, 170910, 501251, 102251,
                348027, 447447, 451859, 151679, 251000, 637019, 307189, 555578),
              start),
  "End" = c(c(154745, 251931, 108824, 460871, 171253, 171027, 501370, 102368,
              348143, 447565, 451975, 151797, 251118, 637136, 307307, 555694),
            end),
  "LenChr" = c(c(228861, 829469, 340914, 1486921, 589812, 299318, 1080440,
                 542723, 449612, 753937, 690901, 1054145, 923535, 791982,
                 1053869, 946846),
               c(219929, 813597, 341580, 1566853, 583092, 271539, 1091538,
                 581049, 440036, 751611, 666862, 1075542, 930506, 777615,
                 1091343, 954457)))


SK1_S288CYuecen <- with(SK1_S288CYuecen,
                        GenomicRanges::GRanges(Chromosome,
                                               IRanges::IRanges(Start, End),
                                               seqlengths=setNames(LenChr,
                                                                   Chromosome)))

SK1_S288CYuecen <- add_genome_name_to_GR(SK1_S288CYuecen, name='SK1_S288C_Yue')



#------------------------------------------------------------------------------#
#                         Intergenic region data                               #
#                         (Conv, div and tandem)                               #
# The intergenic region coordinate file generated using the scripts in:
# '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr2/'
# S288C and old SK1 regions were generated from the files in hwglabr

# 1. Import SK1Yue data
path <- '/Volumes/LabShare/GenomeSequences/hwglabr2'
SK1Yue_intergenic <- read.table(file.path(path, 'SK1Yue_intergenic.txt'),
                                header = TRUE, stringsAsFactors = FALSE)
# 2. Import SK1 data
SK1_intergenic <- read.table(file.path(path, 'SK1_intergenic.txt'),
                             header = TRUE, stringsAsFactors = FALSE)
# 3. Import S288C data
sacCer3_intergenic <- read.table(file.path(path, 'sacCer3_intergenic.txt'),
                                 header = TRUE, stringsAsFactors = FALSE)

#------------------------------------------------------------------------------#
#                           Red1 summits in WT                                 #
#                         (Conv, div and tandem)                               #
# Files generated using the scripts in:
# '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr2/'
# 1. Import SK1Yue data
path <- '/Volumes/LabShare/GenomeSequences/hwglabr2/'
SK1Yue_Red1_summits_file <- file.path(
  path,
  'Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_over20_summits.bed')
SK1Yue_Red1_summits <- rtracklayer::import.bed(SK1Yue_Red1_summits_file)
SK1Yue_Red1_summits <- add_genome_name_to_GR(SK1Yue_Red1_summits, name='SK1Yue')

# 2. Import S288C data
sacCer3_Red1_summits_file <- file.path(
  path,
  'Red1-wildtype-71-34-199-29-Reps-SacCer3-2mis_B3W3_MACS2_over20_summits.bed')
sacCer3_Red1_summits <- rtracklayer::import.bed(sacCer3_Red1_summits_file)
sacCer3_Red1_summits <- add_genome_name_to_GR(sacCer3_Red1_summits,
                                              name='sacCer3')

#------------------------------------------------------------------------------#
#                            Spo11 DSB hotspots                                #
# Files generated using the scripts in:
# '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr2/'
# Creating_spo11_SacCer3_Pan2011_WT1_fixedbegraph.R
# Spo11hotspots_to_SK1Yue.R
# Source of data: nature13120-s2_SacCer2.xls from Thacker 2014 paper.
# S288C data file copied from old package (hwglabr) folder
# Update on 12/2017: New files were made by mapping the Thacker 2014 fastqs to
# the appropriate genomes. New files generated using the script in:
# '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr2/'
# 
path <- '/Volumes/LabShare/GenomeSequences/hwglabr2/'
# 1. Import SK1Yue data
SK1Yue_file <- 'Spo11oligo_WT1_SK1Yue_peaks.bedgraph'
SK1Yue_Spo11_DSBs <- rtracklayer::import.bedGraph(file.path(path, SK1Yue_file))
SK1Yue_Spo11_DSBs <- add_genome_name_to_GR(SK1Yue_Spo11_DSBs, name='SK1Yue')

# 2. Import S288C data
sacCer3_file <- 'Spo11oligo_WT1_Saccer3_peaks.bedgraph'
sacCer3_Spo11_DSBs <- rtracklayer::import.bedGraph(file.path(path,
                                                             sacCer3_file))
sacCer3_Spo11_DSBs <- add_genome_name_to_GR(sacCer3_Spo11_DSBs, name='sacCer3')




#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                           Add all data to package                            #
#                             (as internal data)                               #

# Determine the best compression for the data files
tools::checkRdaFiles('R/') # Suggests 'bzip2'

# Set package directory as working directory
# setwd('/path/to/hwglabr2/')
devtools::use_data(SK1Yuecen, S288CYuecen, sacCer3cen, SK1cen, SK1_S288CYuecen,
                   SK1Yue_intergenic, SK1_intergenic, sacCer3_intergenic,
                   SK1Yue_Red1_summits, sacCer3_Red1_summits,
                   SK1Yue_Spo11_DSBs, sacCer3_Spo11_DSBs,
                   SK1Yue_gff, sacCer3_gff, SK1_gff,
                   internal = TRUE, overwrite = TRUE, compress = "bzip2")
