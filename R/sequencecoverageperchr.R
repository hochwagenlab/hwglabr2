#Function to check Sequence coverage per chromosome
sequencecoverageperchr <- function(dfwithreadcounts,chromosomelength,readlength =350){
  library(reshape2)
  library(tidyverse)
  value = data.frame()
  
  if (missing(dfwithreadcounts)) {
    stop('Please provide a dataframe with readcounts per chromosome ', call.=FALSE)}
  
  if (missing(chromosomelength)) {
    stop('Please provide a length of chromosomes ', call.=FALSE)}
  
  if (missing(readlength)) {
    message("Read length not provided. Using default as 350 bp. ")
    for(i in 1:nrow(dfwithreadcounts)){
      for(j in 1:ncol(dfwithreadcounts)){
        value[i,j] = (dfwithreadcounts[i,j]*readlength)/chromosomelength[i]
      }
    }
  }
  
  else{
    
    for(i in 1:nrow(dfwithreadcounts)){
      for(j in 1:ncol(dfwithreadcounts)){
        value[i,j] = (dfwithreadcounts[i,j]*readlength)/chromosomelength[i]
      }
    }
    colnames(value) = colnames(dfwithreadcounts)
    rownames(value) =rownames(dfwithreadcounts)
    value = melt(rownames_to_column(value))
    
    return(value)
  }
}
