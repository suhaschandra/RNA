library(systemPipeR)
library(GenomicFeatures)


#Read the targets file
targets <- read.delim("targets.txt", comment.char = "#")
#Check if targets file is correct
targets
#Create args
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
moduleload(modules(args))

