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
library(liftOver)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICv2manifest)

# Function to load IDAT files and process methylation data
process_epic_methylation <- function(idat_directory, sample_sheet_path) {
  # Read sample sheet
  sample_sheet <- read.metharray.sheet(sample_sheet_path, pattern = sample_sheet_path)
  
  # Read IDAT files
  raw_data <- read.metharray.exp(base = idat_directory, targets = sample_sheet)
  
  # Preprocess the data (you can choose different normalization methods)
  # Here we'll use functional normalization
  preprocessed_data <- preprocessRaw(raw_data)
  
  # Extract beta values
  beta_values <- getBeta(preprocessed_data)
  
  return(list(
    raw_data = raw_data,
    preprocessed_data = preprocessed_data,
    beta_values = beta_values
  ))
}


extract_mgmt_promoter_betas <- function(beta_values) {
  # mgmt promoter region coordinates (hg38)
  # Adjust these coordinates based on the specific promoter region of interest
  mgmt_promoter_gr <- GRanges(
    seqnames = "chr10",
    ranges = IRanges(start = 131265478, end = 131265604) ## suggested by 10.1186/s13148-021-01044-2
    
  )
  
  # Get annotation data
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  # Convert annotation to GRanges
  annotation_gr <- GRanges(
    seqnames = annotation$chr,
    ranges = IRanges(start = annotation$pos, width = 1),
    strand = "*",
    probe_id = annotation$Name
  )
  
  # Find overlaps between mgmt promoter region and probe locations
  overlaps <- findOverlaps(annotation_gr, mgmt_promoter_gr)
  
  # Extract probe IDs in the mgmt promoter region
  mgmt_promoter_probes <- annotation$Name[queryHits(overlaps)]
  
  # Extract beta values for these probes
  mgmt_promoter_betas <- beta_values[rownames(beta_values) %in% mgmt_promoter_probes, ]
  
  return(list(
    mgmt_promoter_probes = mgmt_promoter_probes,
    mgmt_promoter_betas = mgmt_promoter_betas
  ))
}

# Example usage
# Specify your directories and paths
idat_directory <- "data/D209_D211_D212_D213/2024-386-ILL_METUFB_N=4/METUFB/208789640023"
sample_sheet_path <- "data/D209_D211_D212_D213/2024-386-ILL_METUFB_N=4/samplesheet_2024-386-ILL_METUFB_N=4.csv"

# Process methylation data
methylation_data <- process_epic_methylation(idat_directory, sample_sheet_path)

# Extract mgmt promoter beta values
mgmt_results <- extract_mgmt_promoter_betas(methylation_data$beta_values)

# Print results
print("mgmt Promoter Probes:")
print(mgmt_results$mgmt_promoter_probes)

print("mgmt Promoter Beta Values:")
print(head(mgmt_results$mgmt_promoter_betas))

# Optional: Save results
write.csv(mgmt_results$mgmt_promoter_betas, "mgmt_promoter_beta_values.csv")

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