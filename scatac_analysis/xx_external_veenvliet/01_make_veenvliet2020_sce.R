# Library statements ------------------------------------------------------

library(here)
library(data.table)
library(glue)
library(Matrix)
library(S4Vectors)
library(SingleCellExperiment)

# File specification ------------------------------------------------------

path <- here("data_raw", "external", "veenvliet2020")

barcode_files <- list.files(path, pattern = "barcodes.tsv.gz",
                            recursive = TRUE, full.names = TRUE)

feature_files <- list.files(path, pattern = "features.tsv.gz",
                            recursive = TRUE, full.names = TRUE)

matrix_files <- list.files(path, pattern = "matrix.mtx.gz",
                           recursive = TRUE, full.names = TRUE)

file_struct <- tstrsplit(matrix_files, "/")
samples <- file_struct[[10]]

# Load features -----------------------------------------------------------

features <- lapply(feature_files, fread)


# Load matrices -----------------------------------------------------------

matrices <- lapply(matrix_files, function(f) {
  tmp <- tempfile(fileext = ".mtx")
  cmd <- glue("gunzip -c {f} > {tmp}")
  system(cmd)
  mat <- readMM(tmp)
  mat <- as(mat, "dgCMatrix")
  unlink(tmp)
  return(mat)
})

assay <- do.call(cbind, matrices)

# Load barcodes -----------------------------------------------------------

barcodes <- lapply(barcode_files, fread, header = FALSE)

# Check if same length as matrices
nbarcodes <- vapply(barcodes, nrow, integer(1))
ncols <- vapply(matrices, ncol, integer(1))
all(nbarcodes == ncols) # Should be TRUE

# Convert to coldata
coldat <- rbindlist(barcodes, idcol = "sample")
colDat <- DataFrame(
  barcodes = coldat$V1,
  sample = Rle(samples[coldat$sample])
)


# Features ----------------------------------------------------------------

features <- fread(feature_files[1], header = FALSE)
rowDat <- DataFrame(
  ensembl_id = features$V1,
  gene_symbol = features$V2
)

# Make SCE ----------------------------------------------------------------

sce <- SingleCellExperiment(
  assay = SimpleList(counts = assay),
  colData = colDat,
  rowData = rowDat
)

# Save --------------------------------------------------------------------

saveRDS(sce, here("rds", "XX_veenvliet2020_sce.rds"))
