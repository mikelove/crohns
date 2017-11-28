library(ensembldb)
dbfile <- ensDbFromGtf("Homo_sapiens.GRCh38.90.gtf.gz")
edb <- EnsDb("Homo_sapiens.GRCh38.90.sqlite")
# this is similar to Gencode CHR (except gencode has ~150 additional transcripts in the PAR of Y)
txps <- transcripts(edb, filter=AnnotationFilterList(SeqNameFilter(c(1:22, "X", "Y","MT"))))
library(Biostrings)
# this is missing some of the transcripts in the above
seqs <- readDNAStringSet("Homo_sapiens.GRCh38.cdna.all.fa.gz")
names(seqs) <- sub("(.*?)\\..*","\\1",names(seqs)) # split at '.'
common <- intersect(names(txps), names(seqs)) 
seqs.sub <- seqs[common]
writeXStringSet(seqs.sub, filepath="Homo_sapiens.GRCh38.cdna.std.chroms.fa")
