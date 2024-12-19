# Libraries ---------------------------------------------------------------

library(SingleCellExperiment)
library(GenomicRanges)
library(rtracklayer)
library(here)

# File import -------------------------------------------------------------

out_dir   <- here("data_raw", "bigwig")
exp_file  <- here("rds", "06_LSI_experiment.rds")
read_file <- here("rds", "04_GRanges_joined.rds")

exp   <- readRDS(exp_file)
reads <- readRDS(read_file)

# Filter ------------------------------------------------------------------

read_barcode <- runValue(reads$barcode)
read_batch   <- decode(reads$batch[start(reads$barcode)])
read_alias   <- paste0(read_barcode, "&", read_batch)
reads$alias  <- Rle(read_alias, runLength(reads$barcode))
keep         <- read_alias %in% exp$alias
reads        <- reads[Rle(keep, runLength(reads$barcode))]

reads <- trim(reads) 
reads <- keepStandardChromosomes(reads, "Mus_musculus", "coarse")

# Make bigwigs ------------------------------------------------------------

clusters    <- split(exp$alias, exp$clusters$cluster)
genome_size <- sum(seqlengths(reads))


# All cells ---------------------------------------------------------------

out_file <- file.path(out_dir, "06_All_Timeseries_Cells.bw")

if (!file.exists(out_file)) {
  # Calculate scaling factor to normalise to RPGC
  nbasepairs   <- sum(width(reads))
  scale_factor <- genome_size / nbasepairs
  
  # Compute coverage and scale
  coverage <- coverage(reads)
  coverage <- coverage * scale_factor
  
  # Export
  export.bw(coverage, out_file)
}

# Clusters ----------------------------------------------------------------

for (i in names(clusters)) {
  out_file <- file.path(out_dir, paste0("06_Cluster_", i, ".bw"))
  if (file.exists(out_file)) {
    next
  }
  
  clust_reads  <- reads[reads$alias %in% clusters[[i]]]
  nbasepairs   <- sum(width(clust_reads))
  
  # Calculate scaling factor to normalise to RPGC
  scale_factor <- genome_size / nbasepairs
  
  # Compute coverage and scale
  coverage     <- coverage(clust_reads)
  coverage     <- coverage * scale_factor
  
  # Export
  export.bw(coverage, out_file)
}
