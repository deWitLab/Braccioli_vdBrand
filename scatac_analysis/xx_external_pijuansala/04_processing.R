
# Libraries ---------------------------------------------------------------

library(glue)
library(here)
library(GenomicAlignments)
library(Rsamtools)

# Files -------------------------------------------------------------------

dir <- "/DATA/users/t.vd.brand/projects/mouse_embryo/external_data/sciatac/"

meta_file <- list.files(paste0(dir, "metadata"),
                        pattern = ".xlsx$", full.names = TRUE)
meta <- readxl::read_xlsx(meta_file)

bamfiles  <- list.files(paste0(dir, "bam"),
                        pattern = ".bam$", full.names = TRUE)

blacklist <- rtracklayer::import("/DATA/users/t.vd.brand/projects/dot1l_t_cells/references/blacklist.merge.bed")
blacklist <- as(blacklist, "GNCList")

barcode_out <- paste0(dir, "metadata/barcode_counts.tsv")
out_dir <- paste0(dir, "tabix/")
final_out <- paste0(out_dir, "pijuansala2020.bed")

# Parameters --------------------------------------------------------------

valid_barcodes <- meta$barcode

param <- ScanBamParam(
    flag = scanBamFlag(
        isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE,
        isUnmappedQuery = FALSE, hasUnmappedMate = FALSE
    ),
    mapqFilter = 10L,
    what = "qname"
)

# Functions ---------------------------------------------------------------

shift_fragments <- function(frags, positive = 4, negative = -5) {
    read1 <- tn5shift(GenomicAlignments::first(frags),  positive, negative)
    read2 <- tn5shift(GenomicAlignments::second(frags), positive, negative)
    punion(read1, read2, fill.gap = TRUE, ignore.strand = TRUE)
}

tn5shift <- function(x, positive = 4, negative = -5) {
    x <- resize(granges(x), fix = "start", 1)
    offset <- strand(x) == "+"
    offset <- Rle(
        ifelse(runValue(offset), positive, negative),
        runLength(offset)
    )
    GenomicRanges::shift(x, decode(offset))
}

write_bed <- function(granges, file) {
    granges <- data.table::data.table(
        chrom = as.vector(seqnames(granges)),
        start = start(granges),
        end   = end(granges),
        bc    = granges$barcode
    )
    granges <- granges[, list(duplicates = .N), by = names(granges)]
    data.table::fwrite(
        granges, file,
        append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE
    )
    return(file)
}

read_bam <- function(file) {
    readGAlignmentPairs(file, param = param)
}

clean_barcodes <- function(frags) {
    barcode <- mcols(GenomicAlignments::first(frags))$qname
    subseq(barcode, 1, 22)
}

clean_barcodes_HiSeq2500 <- function(frags) {
    barcode <- mcols(GenomicAlignments::first(frags))$qname
    i1 <- subseq(barcode, 1, 10)
    i2 <- reverseComplement(DNAStringSet(subseq(barcode, 11, 22)))
    as.character(xscat(i1, i2))
}

routine <- function(bamfile, outfile, barcode_out, cleaning,
                    chunksize = 1e6) {
    bamname <- gsub("\\.bam$", "", basename(bamfile))
    bamfile <- BamFile(bamfile, yieldSize = chunksize)
    open(bamfile)
    on.exit(close(bamfile))
    while (length({bamfrag <- read_bam(bamfile)})) {
        barcode <- cleaning(bamfrag)
        bamfrag <- shift_fragments(bamfrag)
        bamfrag$barcode <- barcode
        valid <- barcode %in% valid_barcodes
        valid <- valid & !overlapsAny(bamfrag, blacklist)
        if (mean(valid) < 0.5) {
            message("Less than 50% are valid")
        }
        bamfrag <- bamfrag[valid]
        barcode <- as.data.frame(table(barcode))
        barcode$file <- bamname
        write_bed(bamfrag, outfile)
        data.table::fwrite(
            barcode, barcode_out, append = TRUE,
            quote = FALSE, sep = "\t"
        )
    }
    NULL
}

# Organise data -----------------------------------------------------------

is_HiSeq2500 <- grepl("^HiSeq2500", basename(bamfiles))

names <- gsub("\\.demultiplexed\\.bam$", "", basename(bamfiles))

organise <- list2DF(list(
    bamfiles = bamfiles,
    outfiles = paste0(out_dir, names, ".bed"),
    cleaning = list(clean_barcodes,
                    clean_barcodes_HiSeq2500)[is_HiSeq2500 + 1]
))

# Process files -----------------------------------------------------------

org <- organise[!file.exists(organise$outfiles),]

if (nrow(org) > 0 && !file.exists(final_out)) {

    mapply(routine, bamfile = org$bamfiles, outfile = org$outfiles,
           barcode_out = barcode_out, cleaning = org$cleaning)

    # Summarise barcode counts over all chunks
    barcode_counts <- data.table::fread(barcode_out)
    barcode_counts <- barcode_counts[, .(Freq = sum(Freq)),
                                     by = c("file", "barcode")]
    data.table::fwrite(
        barcode_counts, barcode_out,
        quote = FALSE, sep = "\t"
    )

}

# Combine files -----------------------------------------------------------

if (!file.exists(final_out)) {

    outfiles <- organise$outfiles

    for (file in outfiles) {
        tmp <- data.table::fread(file)
        # Keep only standard chromosomes
        tmp <- tmp[V1 %in% paste0("chr", c(1:19, c("X", "Y")))]
        # Deduplicate
        tmp <- tmp[, .(sum(V5)), by = c("V1", "V2", "V3", "V4")]
        data.table::fwrite(
            tmp, final_out,
            append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE
        )
    }

}

# Make tabix --------------------------------------------------------------

cmd <- glue::glue("sort -k1,1 -k2,2n -o {final_out} {final_out}")
system(cmd)
cmd <- glue::glue("bgzip {final_out}")
system(cmd)
final_out <- paste0(final_out, ".gz")
indexTabix(final_out, seq = 1L, start = 2L, end = 3L)

