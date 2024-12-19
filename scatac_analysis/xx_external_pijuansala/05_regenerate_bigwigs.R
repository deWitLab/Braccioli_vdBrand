# Libraries ---------------------------------------------------------------

library(SingleCellExperiment)
library(GenomicRanges)
library(rtracklayer)
library(here)

# File import -------------------------------------------------------------

out_dir   <- here("data_raw", "external", "pijuansala2020", "bigwigs")
in_dir    <- "/DATA/users/t.vd.brand/projects/mouse_embryo/external_data/sciatac/"
exp_file  <- file.path(in_dir, "metadata", "41556_2020_489_MOESM3_ESM.xlsx")
read_file <- file.path(in_dir, "tabix", "pijuansala2020.bed.gz")

exp   <- readxl::read_excel(exp_file)
reads <- data.table::fread(read_file)

# Filter ------------------------------------------------------------------

seqinfo <- SeqinfoForUCSCGenome("mm10")
seqinfo <- keepStandardChromosomes(seqinfo, "Mus_musculus")
seqinfo <- dropSeqlevels(seqinfo, c("chrY", "chrM"))

reads <- reads[V4 %in% exp$barcode]
reads <- reads[V1 %in% seqnames(seqinfo)]

reads[, cell := match(V4, exp$barcode)]
reads[, cluster := exp$ann[cell]]

reads <- reads[, GRanges(V1, IRanges(V2, V3), cluster = cluster, seqinfo = seqinfo)]
reads <- trim(reads) 

# Make bigwigs ------------------------------------------------------------

clusters    <- setNames(nm = unique(exp$ann))
genome_size <- sum(seqlengths(seqinfo))

for (i in names(clusters)) {
  clust_reads  <- reads[reads$cluster %in% i]
  nbasepairs   <- sum(width(clust_reads))
  
  # Normalise to RPGC
  scale_factor <- genome_size / nbasepairs
  coverage     <- coverage(clust_reads)
  coverage     <- coverage * scale_factor
  
  # Export
  name <- gsub("/| ", "_", i)
  out_file <- file.path(out_dir, paste0(name, ".bw"))
  export.bw(coverage, out_file)
}
