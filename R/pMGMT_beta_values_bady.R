if (F) {
  ## CRAN (updated 21/05/2024)
  install.packages("mgmtstp27_0.8.tar.gz",repos=NULL)
  install.packages("mgmtstp27_0.8.zip",repos=NULL)
  install.packages(c("ade4","MASS"))
  
  ## Bioconductor (https://www.bioconductor.org/install/)
  ## try http:// if https:// URLs are not supported
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("lumi","methylumi","minfi"))
  
  ## test data
  require(mgmtstp27)
  require(minfiData)
  # loading R packages
  # preprocessing of the data
  dat <- preprocessRaw(RGsetEx)
  # computation of M-value
  mvalue <- log2((getMeth(dat)+1)/(getUnmeth(dat)+1))
  mvalue <- as.data.frame(t(mvalue))
  # predictions
  pred1 <- MGMTpredict(mvalue)
  head(pred1)
  # quality control graphics
  par(mfrow=c(2,3))
  MGMTqc.pop(pred1,which.plot=1:3,mfrow=NULL)
  MGMTqc.single(pred1,nsample=1,which.plot=1:3,mfrow=NULL)
} else {
  require(mgmtstp27)
  require(minfiData)
}

## get annotation from: https://support.bioconductor.org/p/9156503/#9156524
library(mgmtstp27)
library(liftOver)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICv2manifest)

# Example usage
# Specify your directories and paths
idat_directory <- "data/D209_D211_D212_D213/2024-386-ILL_METUFB_N=4/METUFB/208789640023"
sample_sheet_path <- "data/D209_D211_D212_D213/2024-386-ILL_METUFB_N=4/samplesheet_2024-386-ILL_METUFB_N=4.csv"

# Read sample sheet
sample_sheet <- read.metharray.sheet(sample_sheet_path, pattern = sample_sheet_path)

# Read IDAT files
raw_data <- read.metharray.exp(base = idat_directory, targets = sample_sheet)

dat <- preprocessRaw(raw_data)
# computation of M-value
mvalue <- log2((getMeth(dat)+1)/(getUnmeth(dat)+1))
mvalue <- as.data.frame(t(mvalue))
# predictions
pred1 <- MGMTpredict(mvalue)
head(pred1)
# quality control graphics
par(mfrow=c(2,3))
MGMTqc.pop(pred1,which.plot=1:3,mfrow=NULL)
MGMTqc.single(pred1,nsample=1,which.plot=1:3,mfrow=NULL)

write.csv(pred1, file.path("results","pred_mgmt.csv"))