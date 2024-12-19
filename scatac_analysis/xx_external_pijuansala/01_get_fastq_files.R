library(here)
library(glue)

dir <- here("fastq/")

sratoolkit <- "/DATA/users/t.vd.brand/tools/sratoolkit.2.9.2-ubuntu64/bin/"

accessions <- readLines(here("fastq", "PRJNA550440_accessions.txt"))

fetch_cmd <- glue("{sratoolkit}prefetch {accessions}")

dump_cmd <- glue("{sratoolkit}fasterq-dump {accessions} -O {dir}")

cmd_gzip <- lapply(accessions, function(acc) {
  acc <- paste0(acc, c("_1.fastq", "_2.fastq"))
  glue("gzip {dir}{acc}")
})

for (i in seq_along(accessions)) {
  system(fetch_cmd[[i]])
  system(dump_cmd[[i]])
  system(cmd_gzip[[i]][[1]])
  system(cmd_gzip[[i]][[2]])
}
