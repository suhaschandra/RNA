# Make sure the right gtf file has been used.
library(systemPipeR)
library(GenomicFeatures)

targets <- read.delim("targets.txt", comment.char = "#")
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")

txdb <- makeTxDbFromGFF(file="data/Macaca_mulatta.Mmul_8.0.1.92.gtf", format="gtf", dataSource="ENSEMBL", organism="Macaca mulatta")
saveDb(txdb, file="./data/Macaca_mulatta8.sqlite")

library(BiocParallel)
txdb <- loadDb("./data/Macaca_mulatta8.sqlite")
eByg <- exonsBy(txdb, by=c("gene"))

# Read the bamfile list
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=TRUE))
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

#Hierarchical Clustering
library(ape)
rpkmDFeByg <- read.delim("./results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[,-19]
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
pdf("results/sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()

#sample outlier analysis using rlog values
library(DESeq2)
countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("results/sample_tree_rlog.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()
                    
#Group-wise PCA using rlog counts
rld <- rlog(dds)
pdf("results/PCA_group_rlog.pdf")
plotPCA(rld)
dev.off()

#Sample-wise PCA using rlog counts
colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
rld <- rlog(dds)
pdf("results/PCA_sample_rlog.pdf")
plotPCA(rld)
dev.off()
                    
#Group-wise PCA based on vsd counts
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_group_vsd.pdf")
plotPCA(vsd)
dev.off()

#Sample-wise PCA based on vsd counts
colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_sample_vsd.pdf")
plotPCA(vsd)
dev.off()

#Groupwise 3D PCA using rlog
library(scatterplot3d)
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
rld <- rlog(dds)
pca <- prcomp(t(assay(rld)), center=TRUE)
pca2 <- pca$x
groups <- targets$Factor
pdf("results/3D_leans.pdf")
scatterplot3d(pca2[,1], pca2[,2], pca2[,3], color = as.numeric(groups), pch=19, type='h')
dev.off()


