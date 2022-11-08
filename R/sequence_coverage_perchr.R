
#' FUNCTION TO CALCULATE SEQUENCING COVERAGE PER CHROMOSOME####

#' INPUTS REQUIRED:
#' dfwithreadcounts (required): Data frame with read counts
#' chromosomelength (required): length of all chromosomes 
#' readlength (optional): Default set to 350 bp


sequence_coverage_perchr <- function(dfwithreadcounts,chromosomelength,readlength =350) {
  t0  <- proc.time()[3]
  library(reshape2)
  library(tidyverse)
  value = data.frame()
  
  if (missing(dfwithreadcounts)) {
    stop('Please provide a dataframe with readcounts per chromosome ', call.=FALSE)}
  
  if (missing(chromosomelength)) {
    stop('Please provide a length of chromosomes ', call.=FALSE)}
 
  else{
    
    for(i in 1:nrow(dfwithreadcounts)){
      for(j in 1:ncol(dfwithreadcounts)){
        value[i,j] = (dfwithreadcounts[i,j]*readlength)/chromosomelength[i]
      }
    }
    colnames(value) = colnames(dfwithreadcounts)
    rownames(value) =rownames(dfwithreadcounts)
    value = melt(rownames_to_column(value))
    
    message('Completed in ', hwglabr2::elapsed_time(t0, proc.time()[3]))
    return(value)
    
  }
  
}
