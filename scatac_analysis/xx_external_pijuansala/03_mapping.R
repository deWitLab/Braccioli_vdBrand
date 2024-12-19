# Library statements ------------------------------------------------------

library(here)
library(glue)

# Set directories and files -----------------------------------------------

path_in  <- here("fastq")
path_out <- here("bam")

files <- list.files(path_in, pattern = ".fastq.gz", full.names = TRUE)

pattern_1 <- "R1.fastq.gz$|R1.repl1.fastq.gz$"
pattern_2 <- "R2.fastq.gz$|R2.repl1.fastq.gz"

reference <- DeWitScripts::refpaths$mouse$mm10$fasta()
read1 <- files[grep(pattern_1, files)]
read2 <- files[grep(pattern_2, files)]

# Check order is correct
correct <- all(gsub(pattern_1, "", read1) == gsub(pattern_2, "", read2))
if (!correct) {
  stop("Order of fastq files is inconsistent.")
}

bamfiles <- gsub(path_in, path_out, read1)
bamfiles <- gsub(pattern_1, "bam", bamfiles)

# Construct commands ------------------------------------------------------

cmd_map <- glue(
  "bwa mem -M -t 20 {reference} {read1} {read2}",
  "samtools view -h -b -q 10",
  "samtools sort -o {bamfiles}",
  .sep = " | "
)

cmd_idx <- glue(
  "samtools index {bamfiles}"
)

# Execute mapping ---------------------------------------------------------

do_these <- !file.exists(bamfiles)

for (i in seq_along(cmd_map)[do_these]) {
  system(cmd_map[[i]])
}

for (i in seq_along(cmd_idx)[do_these]) {
  system(cmd_idx[[i]])
}
