library(here)
sra <- read.delim(here("data","SraRunTable.txt"))
run <- as.character(sra$Run_s)
accession2url <- function(x) {
  prefix <- "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
  dir1 <- paste0("/",substr(x,1,6))
  dir2 <- ifelse(nchar(x) == 9, "",
          ifelse(nchar(x) == 10, paste0("/00",substr(x,10,10)),
          ifelse(nchar(x) == 11, paste0("/0",substr(x,10,11)),
                 paste0("/",substr(x,10,12)))))
  paste0(prefix,dir1,dir2,"/",x)
}
read1 <- file.path(accession2url(run),paste0(run,"_1.fastq.gz"))
read2 <- file.path(accession2url(run),paste0(run,"_2.fastq.gz"))
write(run, here("data","srr"))
write(read1, here("data","read1"))
write(read2, here("data","read2"))

# an example wget call for read1:
#
# > wget -i read1 -nc --random-wait --progress=dot:giga
