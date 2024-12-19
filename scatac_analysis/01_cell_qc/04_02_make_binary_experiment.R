# Note:
# Run as local job with the project directory as working directory

# Libraries ---------------------------------------------------------------

library(GenomicRanges)
library(SummarizedExperiment)
library(Matrix)
library(SingleCellExperiment)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)

# Data --------------------------------------------------------------------

cell_meta <- readRDS("rds/04_cell_metadata.rds")
reads <- readRDS("rds/04_GRanges_joined.rds")
peaks <- readRDS("rds/04_peaks_rough.rds")


# Include only good cells -------------------------------------------------

cell_meta <- cell_meta[cell_meta$QC_stats$inclusion, ]
cell_id <- Rle(paste0(reads$barcode, "&", reads$batch))
reads <- reads[cell_id %in% paste0(as.character(cell_meta$barcode), "&", cell_meta$batch)]
reads <- reads[seqnames(reads) %in% paste0("chr", c(1:19, "X"))]
cell_meta$alias <- paste0(as.character(cell_meta$barcode), "&", cell_meta$batch)
rownames(cell_meta) <- cell_meta$alias

# Make coverage track -----------------------------------------------------

#cover <- coverage(reads)
#export(cover, "../../processed_data/sciatac/bigwig/04_good_cells.bw", 
#       format = "BigWig")


# Calculate CG content ----------------------------------------------------

cell_id <- Rle(paste0(reads$barcode, "&", reads$batch))
cells <- split(reads, cell_id)
GC_content <- bplapply(cells, function(cell) {
  seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, cell)
  abcfreq <- colSums(alphabetFrequency(seq))
  sum(abcfreq[2:3]) / sum(abcfreq[1:4])
}, BPPARAM = MulticoreParam(workers = 10))
GC_content <- unlist(GC_content)
cell_meta$QC_stats$GC_content <- GC_content[cell_meta$alias]

# Make binary matrix ------------------------------------------------------

# Split reads to barcodes and reduce ranges

# cells <- reduce(cells)

# Make matrix
nr <- NROW(peaks)
peaks <- as(peaks, "GNCList")
ncells <- length(cells)

# Chunking over cells, since vapply returns a dense matrix (memory problem)
chunksize <- 1000
chunks <- ceiling(ncells / chunksize)

i <- split(seq_len(ncells), rep(seq_len(chunks), each = chunksize)[seq_len(ncells)])
mats <- lapply(i, function(ii) {
  matrix <- vapply(cells[ii], function(cell) {
    overlapsAny(peaks, cell)
  }, logical(nr))
  matrix <- as(matrix, "dgCMatrix")
  matrix
})
mats <- do.call(cbind, mats)
colnames(mats) <- names(cells)

# Synchronise order of cell metadata
cell_meta <- cell_meta[match(colnames(mats), cell_meta$alias),]

# Make experiment
full_exp <- SingleCellExperiment(
  assays = list(binary = mats),
  rowRanges = as(peaks, "GRanges"),
  colData = cell_meta
)

# Capture metadata --------------------------------------------------------

generator_script <- as.character(sys.call(1))[2]
generator_script <- if (!exists("generator_script")) "manual" else generator_script
system_info <- Sys.info()
session_info <- utils::sessionInfo()

metadata <- list(
  intermediate_file = list(normalizePath("rds/04_peaks_rough.rds"),
                           normalizePath("rds/04_GRanges_joined.rds"),
                           normalizePath("rds/04_cell_metadata.rds")),
  generator = generator_script,
  system = system_info,
  session = session_info
)

# Write files -------------------------------------------------------------

metadata(full_exp) <- metadata
saveRDS(full_exp, file = "rds/04_binary_experiment.rds")
