library(tidyverse)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
#library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(GEOquery)

## load data
options(timeout = max(30000000000000000000, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

accessions <- c("GSE164994","GSE101638","GSE139652")
map(accessions[2], function(i) { #[159:length(accessions)] 85 ans 93,94,95,96 did not work
  #get raw data - idats, processed beta matrix, etc.
  #dir.create(file.path("data",i))
  #
  tryCatch({
  getGEOSuppFiles(i, baseDir = file.path("data"))
  gse <- getGEO(i, GSEMatrix = TRUE, destdir = file.path("data",i))
  dir.create(file.path("data",i,"idat"))
  untar(file.path("data",i, paste0(i,"_RAW.tar")), exdir = file.path("data",i,"idat"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
