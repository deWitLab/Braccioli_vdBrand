# Note:
# Run as local job with the project directory as working directory

# Libraries ---------------------------------------------------------------

library(ShortRead)
library(Biostrings)

# Set directories and files -----------------------------------------------

path_in <- "../../raw_data/sciatac/fastq/"
path_out <- normalizePath("../../processed_data/sciatac/fastq_namepaste/")

files <- list.files(path_in, pattern = "fastq.gz", full.names = TRUE)

read1 <- normalizePath(files[grepl("_R1_001.fastq.gz", files)])
read2 <- normalizePath(files[grepl("_R2_001.fastq.gz", files)])

# Handle read 1 -----------------------------------------------------------

stream1 <- FastqStreamer(read1)

while (length(fq1 <- yield(stream1))) {
  fq1_id <- id(fq1)
  nchar <- width(fq1_id)
  # Grab barcode
  barcode <- subseq(fq1_id, nchar - 32, nchar)
  # Filter out the '+' symbol
  barcode <- xscat(subseq(barcode, 1, 16), subseq(barcode, 18, 33))
  # Prepend SAM specification barcode identifiers
  barcode <- xscat("BC:Z:", barcode)
  
  # Grab readname
  new_id <- xscat(subseq(fq1_id, 1, nchar - 39), barcode)
  fq1@id <- new_id
  writeFastq(fq1, file = paste0(path_out, "/Undetermined_R1_001.fastq.gz"),
             mode = "a", compress = TRUE)
}

close(stream1)

# Handle read 2 -----------------------------------------------------------

stream2 <- FastqStreamer(read2)

while (length(fq2 <- yield(stream2))) {
  fq2_id <- id(fq2)
  nchar <- width(fq2_id)
  # Grab barcode
  barcode <- subseq(fq2_id, nchar - 32, nchar)
  # Filter out the '+' symbol
  barcode <- xscat(subseq(barcode, 1, 16), subseq(barcode, 18, 33))
  # Prepend SAM specification barcode identifiers
  barcode <- xscat("BC:Z:", barcode)
  
  # Grab readname
  new_id <- xscat(subseq(fq2_id, 1, nchar - 39), barcode)
  fq2@id <- new_id
  writeFastq(fq2, file = paste0(path_out, "/Undetermined_R2_001.fastq.gz"),
             mode = "a", compress = TRUE)
}

close(stream2)
