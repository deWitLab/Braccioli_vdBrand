# Note:
# Run as local job with the project directory as working directory

# Libraries ---------------------------------------------------------------

library(ShortRead)
library(Biostrings)
library(here)

# Set directories and files -----------------------------------------------

in_dir <- here("fastq")
out_dir <- here("fastq", "processed")

files <- list.files(in_dir, pattern = "fastq.gz", full.names = TRUE)
outfiles <- gsub(in_dir, out_dir, files)

# Function ----------------------------------------------------------------

# Essentially we want to cut out the cellular barcode bit of the read header
# and format it such that bwa mem can recognise it.

process_fastq <- function(file, outfile, sep = ":", place = 1) {
  stream <- FastqStreamer(file)
  on.exit(close(stream))
  
  while (length(reads <- yield(stream))) {
    id <- id(reads)
    id <- as.character(id)
    id <- data.table::tstrsplit(id, sep)
    
    barcode <- id[[place]]
    id <- do.call(paste, c(id[-place], sep = sep))
    id <- paste0(id, " BC:Z:", barcode)
    id <- as(id, "BStringSet")
    
    reads@id <- id
    writeFastq(reads, file = outfile, mode = "a", compress = TRUE)
  }
}

# Process -----------------------------------------------------------------

for (i in seq_along(files)) {
  if (file.exists(outfiles[[i]])) {
    next()
  }
  process_fastq(files[[i]], outfiles[[i]])
}
