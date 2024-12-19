
library(here)
library(data.table)
library(GenomicAlignments)
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)

primers <- here("rds", "primers_timecourse.rds") 
primers <- readRDS(primers)

theobarcode <- split(primers$sequence, primers$type)
theobarcode <- list(
  "P7" = as.character(theobarcode$P7),
  "P5" = as.character(subseq(reverseComplement(DNAStringSet(theobarcode$P5)), 2, 8)),
  "i7" = as.character(subseq(theobarcode$i7, 1, 7)),
  "i5" = as.character(reverseComplement(subseq(DNAStringSet(theobarcode$i5), 1, 7)))
)

bamfiles <- here("data_raw", "bam", "kos")
bamfiles <- list.files(bamfiles, pattern = ".bam$", full.names = TRUE)
bamfiles <- BamFile(bamfiles, yieldSize = 1e6)

bed_temp <- tempfile(fileext = ".bed")
barcode_temp <- tempfile(fileext = ".tsv")
seqinf   <- SeqinfoForUCSCGenome("mm10")

param <- ScanBamParam(tag = "BC")

open(bamfiles)
while (length({bamfrag <- readGAlignmentPairs(bamfiles, param = param)})) {
  bamfrag <- bamfrag[bamfrag@isProperPair]
  barcode <- mcols(GenomicAlignments::first(bamfrag))$BC
  
  pieces <- list(
    "P7" = substr(barcode, 2,  9),
    "P5" = substr(barcode, 12, 18),
    "i7" = substr(barcode, 21, 27),
    "i5" = substr(barcode, 29, 40)
  )
  
  matches <- Map(`%in%`, x = pieces, table = theobarcode)
  matches <- Reduce(`+`, matches)
  valid   <- matches == 4
  
  bamfrag <- bamfrag[valid]
  
  mate1 <- granges(GenomicAlignments::first(bamfrag, real.strand = TRUE))
  mate2 <- granges(GenomicAlignments::last(bamfrag, real.strand = TRUE))
  mate1 <- IRanges::shift(mate1, ifelse(strand(mate1) == "+", 4, -5))
  mate2 <- IRanges::shift(mate2, ifelse(strand(mate2) == "+", 4, -5))
  
  frags <- punion(mate1, mate2, fill.gap = TRUE, ignore.strand = TRUE)
  frags <- data.table(
    chrom = as.vector(seqnames(frags)),
    start = start(frags),
    end   = end(frags),
    strand = as.vector(strand(frags)),
    bc = mcols(GenomicAlignments::first(bamfrag))$BC
  )
  frags <- frags[, .(duplicates = .N), by = names(frags)]
  
  fwrite(frags, bed_temp, append = TRUE, quote = FALSE, sep = "\t")
  fwrite(as.data.table(table(barcode)),
         barcode_temp, append = TRUE, quote = FALSE, sep = "\t")
  
}
close(bamfiles)

# Format fragments --------------------------------------------------------

allreads <- fread(bed_temp)
allreads <- allreads[, list(duplicates = sum(duplicates)), 
                     by = eval(head(names(allreads), -1))]

allreads <- allreads[, GRanges(chrom, IRanges(start, end), strand,
                               barcode = bc, duplicates = duplicates,
                               seqinfo = seqinf)]

saveRDS(allreads, here("rds", "06_GRanges_KO.rds"))

# Format barcode stats ----------------------------------------------------

bcs <- fread(barcode_temp)
bcs <- bcs[, sum(N), by = c("barcode")]
saveRDS(bcs, here("rds", "06_barcode_counts_KO.rds"))
