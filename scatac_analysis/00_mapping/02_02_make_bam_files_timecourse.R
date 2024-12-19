# Note:
# Run as local job with the project directory as working directory

# Set directories and files -----------------------------------------------

path_in <- normalizePath("../../processed_data/sciatac/fastq_namepaste/timecourse/")
path_out <- normalizePath("../../processed_data/sciatac/bam/timecourse/")

files <- list.files(path_in, pattern = "fastq.gz", full.names = TRUE)

reference <- normalizePath("/DATA/references/mouse/mm10/mm10.fa")
read1 <- normalizePath(files[grepl("_R1_001.fastq.gz", files)])
read2 <- normalizePath(files[grepl("_R2_001.fastq.gz", files)])
bamfile <- paste0(path_out, "/Undetermined_0_sorted.bam")

# Call mapper -------------------------------------------------------------

command <- paste0("bwa mem -MC -t 18 ",
                  reference, " ", read1, " ", read2,
                  " | samtools view -h -b -q 10 | samtools sort -o ",
                  bamfile)

edit_tags <- paste0("samtools view -h ", bamfile, " | ",
                    "awk -F'BC:Z:' '{ if($0 ~ ", '"^@"',
                    ") {print $0} else {print $1 ",
                    '"BC:Z: "',
                    "substr($2, 1, 23) substr($2, 25, 7)} }' | ", "samtools view -b -o ",
                    path_out, "timecourse_tagedit.bam")
edit_tags <- gsub("\\\\", "", edit_tags, perl = TRUE)

system(command)
system(paste0("samtools index ", bamfile))
