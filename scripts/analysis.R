library(here)
library(readr)
coldata <- read_delim(here("data","SraRunTable.txt"), delim="\t")
files <- here("data","quants",coldata$Run_s,"quant.sf.gz")
names(files) <- coldata$Run_s
all(file.exists(files))

suppressPackageStartupMessages(library(GenomicFeatures))
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
# takes 100 seconds to make TxDb
# gtf <- makeTxDbFromGFF("gencode.v27.annotation.gtf.gz")
# saveDb(gtf, file="gencode.v27.sqlite")
gtf <- loadDb("gencode.v27.sqlite")
columns(gtf)
tx2gene <- select(gtf, keys(gtf, "TXNAME"), "GENEID", "TXNAME")

library(tximport)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# this is always necessary...
coldata$disease_stage_s
coldata$condition <- coldata$disease_stage_s
coldata$condition[is.na(coldata$condition)] <- "normal"
coldata$condition <- sub("/.*","",coldata$condition)
coldata$condition <- factor(coldata$condition, c("normal","B1","B2","B3"))

suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)

# remove failing samples
table(dds$condition)
qcFail <- c("SRR1813874","SRR1813893")
qcFail %in% colnames(dds)
dds <- dds[,!colnames(dds) %in% qcFail]
table(dds$condition)

# remove other samples
otherRemove <- c("SRR1813895","SRR1813889")
otherRemove %in% colnames(dds)
dds <- dds[,!colnames(dds) %in% otherRemove]
table(dds$condition)


vsd <- vst(dds)
plotPCA(vsd)
plotPCA(vsd, "Sex_s")
plotPCA(vsd, "age_s")
plotPCA(vsd, "MBases_l")

## suppressPackageStartupMessages(library("sva"))
## dds <- estimateSizeFactors(dds)
## dat  <- counts(dds, normalized=TRUE)
## idx  <- rowMeans(dat) > 1
## dat  <- dat[idx, ]
## mod  <- model.matrix(~ condition, colData(dds))
## mod0 <- model.matrix(~   1, colData(dds))
## fit <- svaseq(dat, mod, mod0, n.sv = 2)
## dds$sv1 <- fit$sv[,1]
## dds$sv2 <- fit$sv[,2]
## vsd$sv1 <- fit$sv[,1]
## vsd$sv2 <- fit$sv[,2]
## plotPCA(vsd, "sv1")
## plotPCA(vsd, "sv2")

keep <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 5
table(keep)
dds <- dds[keep,]

# takes 1.5 min
design(dds) <- ~condition
system.time({ dds <- DESeq(dds, test="LRT", reduced=~1, minReplicatesForReplace=Inf) })
res <- results(dds, alpha=.05)
summary(res)

par(mfrow=c(3,3))
for (i in 1:9) 
  plotCounts(dds, gene=order(res$pvalue)[i], transform=FALSE)

rownames(res)[order(res$pvalue)[1:9]]

library(org.Hs.eg.db)
unname(
  mapIds(org.Hs.eg.db,
       sub("\\..*","",rownames(res)[order(res$pvalue)[1:50]]),
       "SYMBOL",
       "ENSEMBL")
)

# Ben # normal samples vs 
