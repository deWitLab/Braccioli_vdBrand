# Libraries ---------------------------------------------------------------

library(GenomicRanges)
library(SummarizedExperiment)
library(Matrix)
library(SingleCellExperiment)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table, exclude = "shift")
library(here)

# Load data ---------------------------------------------------------------

cell_meta <- readRDS(here("rds", "05_cell_metadata_KOs.rds"))
reads     <- readRDS(here("rds", "06_GRanges_KO.rds"))
peaks     <- readRDS(here("rds", "04_peaks_rough.rds"))

# Include only good cells -------------------------------------------------

cell_meta <- cell_meta[cell_meta$inclusion, ]
reads <- reads[reads$barcode %in% cell_meta$barcode]
reads <- reads[seqnames(reads) %in% paste0("chr", c(1:19, "X"))]
cell_meta$alias <- paste0(cell_meta$barcode, "&", cell_meta$batch)
rownames(cell_meta) <- cell_meta$alias

# Calculate CG content ----------------------------------------------------

seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, reads)
seq <- alphabetFrequency(seq)
seq <- cbind(rowSums(seq[, c("C", "G")]), rowSums(seq[, c("A", "T")]))
seq <- split.data.frame(seq, reads$barcode)
seq <- vapply(seq, colSums, numeric(2))
seq <- t(seq)
seq <- setNames(seq[, 1] / rowSums(seq), rownames(seq))

cell_meta$QC_stats$GC_content <- seq[cell_meta$barcode]

# Make matrix -------------------------------------------------------------

fo <- findOverlaps(reads, peaks)
fo <- data.table(barcode = reads$barcode[from(fo)],
                 peak    = to(fo))
fo <- fo[, .N, by = c("barcode", "peak")]
fo[, match := match(barcode, cell_meta$barcode)]
mdim <- c(length(peaks), nrow(cell_meta))

matrix <- fo[, sparseMatrix(peak, match, x = rep(1, length(peak)), dims = mdim)]

full_exp <- SingleCellExperiment(
  assays = list(binary = matrix),
  rowRanges = peaks,
  colData = cell_meta
)

saveRDS(full_exp, here("rds", "05_binary_experiment_KOs.rds"))
