pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE,repos="http://bioconductor.org/biocLite.R")
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

pkgTest("seqinr")
args <- commandArgs(TRUE)
script<-args[1]
motif<-args[2]
mod_shift<-as.numeric(args[3])
ref<-args[4]

setwd("./")
source(script)
findMotifReadSpace(motif=motif,shift=mod_shift,genomeName=ref)
