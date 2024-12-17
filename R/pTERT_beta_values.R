# Required libraries
library(minfi)
library(GenomicRanges)

# Function to load IDAT files and process methylation data
process_epic_methylation <- function(idat_directory, sample_sheet_path) {
  # Read sample sheet
  sample_sheet <- read.metharray.sheet(idat_directory, pattern = sample_sheet_path)
  
  # Read IDAT files
  raw_data <- read.metharray.exp(base = idat_directory, targets = sample_sheet)
  
  # Preprocess the data (you can choose different normalization methods)
  # Here we'll use functional normalization
  preprocessed_data <- preprocessFunNorm(raw_data)
  
  # Extract beta values
  beta_values <- getBeta(preprocessed_data)
  
  return(list(
    raw_data = raw_data,
    preprocessed_data = preprocessed_data,
    beta_values = beta_values
  ))
}

# Function to extract beta values for TERT promoter region
tert_atg <- 1295104
prom_start <- tert_atg+2300 ## upstream based on https://link.springer.com/article/10.1007/s00432-021-03536-3

extract_tert_promoter_betas <- function(beta_values) {
  # TERT promoter region coordinates (hg38)
  # Adjust these coordinates based on the specific promoter region of interest
  tert_promoter_gr <- GRanges(
    seqnames = "chr5",
    #ranges = IRanges(start = 1295005, end = 1296876) ## suggested by claude ai
    ranges = IRanges(start = tert_atg, end = prom_start)
  )
  
  # Get annotation data
  annotation <- getAnnotation(Epic850k.hg38)
  
  # Convert annotation to GRanges
  annotation_gr <- GRanges(
    seqnames = annotation$chr,
    ranges = IRanges(start = annotation$pos, width = 1),
    strand = "*",
    probe_id = annotation$Name
  )
  
  # Find overlaps between TERT promoter region and probe locations
  overlaps <- findOverlaps(annotation_gr, tert_promoter_gr)
  
  # Extract probe IDs in the TERT promoter region
  tert_promoter_probes <- annotation$Name[queryHits(overlaps)]
  
  # Extract beta values for these probes
  tert_promoter_betas <- beta_values[rownames(beta_values) %in% tert_promoter_probes, ]
  
  return(list(
    tert_promoter_probes = tert_promoter_probes,
    tert_promoter_betas = tert_promoter_betas
  ))
}

# Example usage
# Specify your directories and paths
idat_directory <- "/path/to/your/idat/files/"
sample_sheet_path <- "sample_sheet.csv"

# Process methylation data
methylation_data <- process_epic_methylation(idat_directory, sample_sheet_path)

# Extract TERT promoter beta values
tert_results <- extract_tert_promoter_betas(methylation_data$beta_values)

# Print results
print("TERT Promoter Probes:")
print(tert_results$tert_promoter_probes)

print("TERT Promoter Beta Values:")
print(head(tert_results$tert_promoter_betas))

# Optional: Save results
write.csv(tert_results$tert_promoter_betas, "tert_promoter_beta_values.csv")