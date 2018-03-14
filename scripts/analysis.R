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
lvls <- c("normal","B1","B2","B3")
coldata$condition <- factor(coldata$condition, lvls)

suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)

# remove failing samples
table(dds$condition)
qcFail <- c("SRR1813874","SRR1813893")
qcFail %in% colnames(dds)
dds <- dds[,!colnames(dds) %in% qcFail]
table(dds$condition)

vsd <- vst(dds)
plotPCA(vsd)
plotPCA(vsd, "Sex_s")
plotPCA(vsd, "age_s")
plotPCA(vsd, "MBases_l")

# add new metadata
gsm <- read_delim(here("data","GSM_table.tsv"),delim="\t",col_names=c("GSM","title"))
samp <- read_delim(here("data","sample_table.tsv"),delim="\t")
gsm$sample <- sub(".*\\[RNA-Seq, (.*?)_.*\\]","\\1",gsm$title)
samp$GSM <- gsm$GSM[match(samp$sample, gsm$sample)]
for (column in c("tissue", "montreal", "batch")) {
  colData(dds)[[column]] <- samp[[column]][match(dds$Sample_Name_s, samp$GSM)]
}
dds$montreal <- factor(dds$montreal, lvls)
dds <- dds[,dds$tissue == "colon"]
table(dds$condition, dds$montreal)

vsd <- vst(dds)
plotPCA(vsd)
plotPCA(vsd, "montreal")
plotPCA(vsd, "batch")

design(dds) <- ~montreal

dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 5
table(keep)
dds <- dds[keep,]

# takes <1 min
system.time({ dds <- DESeq(dds, test="LRT", reduced=~1, minReplicatesForReplace=Inf) })
res <- results(dds, alpha=.05)
summary(res)

par(mfrow=c(3,3))
for (i in 1:9) 
  plotCounts(dds, gene=order(res$pvalue)[i], transform=FALSE)
par(mfrow=c(1,1))

rownames(res)[order(res$pvalue)[1:9]]

library(org.Hs.eg.db)
unname(
  mapIds(org.Hs.eg.db,
       sub("\\..*","",rownames(res)[order(res$pvalue)[1:50]]),
       "SYMBOL",
       "ENSEMBL")
)

library(goseq)
padj <- ifelse(is.na(res$padj), 1, res$padj)
genes <- as.integer(padj < .05)
strp <- function(x) sub("(.*)\\..*","\\1",x)
names(genes) <- strp(rownames(res))
table(genes)
table(duplicated(names(genes)))
table(genes[duplicated(names(genes))])
genes <- genes[!duplicated(names(genes))]

supportedOrganisms()[supportedOrganisms()$Genome=="hg38",]

# it's ok! we can do this
ebg <- exonsBy(gtf, by="gene")
rangePlotter <- function(x) {
  plot(c(start(x),end(x)),rep(seq_along(x),2),cex=.5,xlab="bp",ylab="")
  segments(start(x),seq_along(x),end(x),seq_along(x),lwd=3)
}
rangePlotter(ebg[[1]])
ebg.red <- reduce(ebg)
rangePlotter(ebg.red[[1]])
lengthData <- sum(width(ebg))
names(lengthData) <- strp(names(lengthData))
lengthData <- lengthData[names(genes)]

pwf <- nullp(genes, genome="hg39", id="ensGene", bias.data=lengthData)
go <- goseq(pwf, "hg38", "ensGene", test.cats="GO:BP")

head(subset(go, numInCat < 200),20)

