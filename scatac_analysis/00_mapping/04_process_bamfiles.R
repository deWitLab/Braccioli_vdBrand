library(data.table)
library(GenomicAlignments)
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)

primers <- "rds/primers_timecourse.rds"
primers <- readRDS(primers)

bamfiles <- list.files(here::here("data_raw", "bam"), recursive = TRUE,
                       pattern = "tagedit.bam$", full.names = TRUE)
batches <- tstrsplit(bamfiles, "/")[[9]]

# Determine valid barcodes ------------------------------------------------

## Note that we need to shorten and reverse-complement the i5 and P5 primers,
## otherwise we won't find any matches within the read names
valid_barcodes <- as(ifelse(
  primers$type %in% c("i7", "P7"),
  primers$sequence,
  reverseComplement(subseq(primers$sequence, 1, 7))
), "DNAStringSet")

valid_barcodes <- split(valid_barcodes, primers$type)
valid_barcodes <- expand.grid("i7" = valid_barcodes$i7, 
                              "P7" = valid_barcodes$P7,
                              "i5" = valid_barcodes$i5, 
                              "P5" = valid_barcodes$P5)
valid_barcodes <- do.call(xscat, valid_barcodes)

# Initialise --------------------------------------------------------------

bed_temp <- tempfile(fileext = ".bed")
barcode_temp <- tempfile(fileext = ".tsv")
seqinf <- lapply(bamfiles, function(x){seqinfo(BamFile(x))})
seqinf <- do.call(merge, seqinf)
genome(seqinf) <- "mm10"

param <- ScanBamParam(tag = "BC")
i <- 1

# Loop through bam files --------------------------------------------------

for(j in seq_along(bamfiles)) {
  bamfile <- BamFile(bamfiles[[j]], yieldSize = 1e6)
  open(bamfile)
  while(length({bamfrag <- readGAlignmentPairs(bamfile, param = param)})) {
    # Filter inproper pairs
    bamfrag <- bamfrag[bamfrag@isProperPair]
    barcode <- mcols(GenomicAlignments::first(bamfrag))$BC
    if ({len <- nchar(barcode[1])} == 31) {
      barcode <- subseq(barcode, 2, len)
    }
    valid <- barcode %in% valid_barcodes
    bamfrag <- bamfrag[valid]
    
    # Apply Tn5 shift
    mate1 <- granges(GenomicAlignments::first(bamfrag, real.strand = TRUE))
    mate2 <- granges(GenomicAlignments::last(bamfrag, real.strand = TRUE))
    mate1 <- IRanges::shift(mate1, ifelse(strand(mate1) == "+", 4, -5))
    mate2 <- IRanges::shift(mate2, ifelse(strand(mate2) == "+", 4, -5))
    
    # Deduplicate
    frags <- punion(mate1, mate2, fill.gap = TRUE, ignore.strand = TRUE)
    frags <- data.table(
      chrom = as.vector(seqnames(frags)),
      start = start(frags),
      end = end(frags),
      strand = as.vector(strand(frags)),
      bc = barcode[valid],
      batch = batches[j]
    )
    frags <- frags[, list(duplicates = .N), by = names(frags)]
    
    # Save info
    fwrite(frags, bed_temp, append = TRUE, quote = FALSE, sep = "\t")
    fwrite(as.data.table(table(barcode))[, batch := batches[j]], 
           barcode_temp, append = TRUE,
           quote = FALSE, sep = "\t")
  }
  close(bamfile)
}



# Format fragments --------------------------------------------------------

allreads <- fread(bed_temp)
allreads <- allreads[, list(duplicates = sum(duplicates)), 
                     by = eval(head(names(allreads), -1))]
setkeyv(allreads, c("batch", "bc", "chrom"))

allreads <- allreads[, GRanges(chrom, IRanges(start, end), strand,
                               barcode = bc, duplicates = duplicates,
                               batch = batch,
                               seqinfo = seqinf)]
allreads$barcode <- Rle(allreads$barcode)
allreads$batch <- Rle(allreads$batch)

saveRDS(allreads, "rds/04_GRanges_joined.rds")

# Format barcode stats ----------------------------------------------------

bcs <- fread(barcode_temp)
bcs <- bcs[, sum(N), by = c("barcode", "batch")]
saveRDS(bcs, "rds/04_barcode_counts_joined.rds")
