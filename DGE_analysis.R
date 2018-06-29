#Run egdeR following read counting and quality control

library(systemPipeR)
library(GenomicFeatures)
library(edgeR)

#Read the targets file
targets <- read.delim("targets.txt", comment.char = "#")
countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")

#For monkeys
desc <- read.delim("/bigdata/messaoudilab/ssure003/Reference/Rhesus_annotations.xls", row.names=1)
#For humans follow the following
desc <- read.delim("/bigdata/messaoudilab/ssure003/Reference/Homo_sapiens/Human_Annotation.txt", row.names=1)

edgeDF <- cbind(edgeDF, desc[rownames(edgeDF),])
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)

edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
pdf("results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=5))
dev.off()
