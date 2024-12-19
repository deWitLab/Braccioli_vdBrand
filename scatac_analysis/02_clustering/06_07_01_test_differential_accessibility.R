# Libraries ---------------------------------------------------------------

library(SingleCellExperiment)
library(here)
library(Matrix)
library(BiocParallel)

# Load data ---------------------------------------------------------------

exp_file <- here("rds", "06_LSI_experiment.rds")
exp <- readRDS(exp_file)

# Settings ----------------------------------------------------------------

out_file <- here("rds", "06_peak_statistics.rds")

# Proportion of cells to sample from every cluster as reference level
ref_prop  <- 0.2

# Minimum proportion of cells containing a peak
cell_prop <- 0.02

# Filter ------------------------------------------------------------------

peaks <- rowRanges(exp)
peaks$index <- seq_along(peaks)

mat <- assay(exp, "counts")

keep_cells <- exp$batch != "KOs"

row_means  <- rowMeans(mat[, keep_cells])
keep_peaks <- row_means >= cell_prop

exp   <- exp[keep_peaks, keep_cells]
peaks <- peaks[keep_peaks] 

# Setup data --------------------------------------------------------------

## Accessibility ----------------------------------------------------------

# For every peak, keep indices of cells that accessible there
# This is mostly to avoid going to memory-hungry dense data
mat  <- assay(exp, "counts")
mat  <- as(mat, "dgTMatrix")
data <- split(mat@j + 1, mat@i + 1)

## Covariates -------------------------------------------------------------

gc_bias   <- unname(exp$QC_stats$GC_content)
cell_size <- unname(colMeans(mat))

# Set largest batch as reference
batch <- factor(unname(exp$batch))
batch <- relevel(batch, ref = tail(names(sort(table(batch))), 1))

## Reference cells --------------------------------------------------------

set.seed(20220415)
cluster   <- as.character(exp$clusters$cluster)

# Sample cells from every cluster for a reference level
ref <- split(seq_along(cluster), cluster)
ref <- lapply(ref, function(x) {sample(x, size = round(length(x) * ref_prop))})
ref <- unlist(ref)

# Set reference level
cluster[ref] <- "Ref"
cluster      <- factor(cluster, levels = sort(unique(cluster)))
cluster      <- relevel(cluster, ref = "Ref")

## Template ---------------------------------------------------------------

# We make a template dataframe where we just have to swap out `y`
# for every peak, but can retain the rest of data
template <- data.frame(
  y       = 0,
  size    = cell_size,
  gc_bias = gc_bias,
  batch   = batch,
  cluster = cluster
)

# Tests -------------------------------------------------------------------

## Setup models -----------------------------------------------------------

RhpcBLASctl::blas_set_num_threads(1)

formulae <- list(
  full    = y ~ 1 + size + gc_bias + batch + cluster,
  reduced = y ~ 1 + size + gc_bias + batch 
)

models      <- lapply(formulae, model.matrix, data = template)
coefnames   <- paste0("mcluster", sort(unique(cluster)))[-1]

# Names for output
prettynames <- gsub("cluster", "cluster_", gsub("^m", "", coefnames))
est_names   <- paste0(prettynames, "_", "estimate")
pval_names  <- paste0(prettynames, "_", "pvalue")
ll_names    <- paste0("loglik_", names(formulae))
out_names   <- c(est_names, pval_names, "intercept", ll_names)

## Run tests --------------------------------------------------------------

tests <- bplapply(data, function(i) {
  
  df <- template
  
  # Set `y` to 1 where we have cell indices
  df$y[i] <- 1
  
  # Fit models
  fits <- lapply(models, function(m) {
    glm(y ~ m, data = df, family = "binomial")
  })
  
  logliks <- vapply(fits, logLik, numeric(1))
  coefs   <- summary(fits$full)$coefficients
  
  out <- c(as.vector(coefs[coefnames, c(1, 4)]), coefs[1, 1], logliks)
  
  setNames(out, out_names)
}, BPPARAM = MulticoreParam(10))

## Format output ----------------------------------------------------------

tests <- do.call(rbind, tests)
tests <- as(tests, "DataFrame")
tests$ranges <- peaks

## Likelihood ratio test ---------------------------------------------------

df <- vapply(models, ncol, integer(1))
df <- unname(abs(diff(df)))

tests$lr_stat <- -2 * (tests$loglik_reduced - tests$loglik_full)
tests$lr_pval <- pchisq(tests$lr_stat, df = df, lower.tail = FALSE)


# Save results ------------------------------------------------------------

saveRDS(tests, out_file)

