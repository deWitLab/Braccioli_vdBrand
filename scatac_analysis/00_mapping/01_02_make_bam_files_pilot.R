# Note:
# Run as local job with the project directory as working directory

# Set directories and files -----------------------------------------------

path_in <- normalizePath("../../processed_data/sciatac/fastq_namepaste/")
path_out <- normalizePath("../../processed_data/sciatac/bam/")

files <- list.files(path_in, pattern = "fastq.gz", full.names = TRUE)

reference <- normalizePath("/DATA/references/mouse/mm10/mm10.fa")
read1 <- normalizePath(files[grepl("_R1_001.fastq.gz", files)])
read2 <- normalizePath(files[grepl("_R2_001.fastq.gz", files)])
bamfile <- paste0(path_out, "/Undetermined_0_sorted.bam")

# Call mapper -------------------------------------------------------------

command <- paste0("bwa mem -MC -t 20 ",
                  reference, " ", read1, " ", read2,
                  " | samtools view -h -b -q 10 | samtools sort -o ",
                  bamfile)

system(command)
system(paste0("samtools index ", bamfile))
