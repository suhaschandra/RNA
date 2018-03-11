#The following code has been modified from Dr Girke's systemPipeR pipeline.

library(systemPipeR)
library(GenomicFeatures)
#Read the targets file
targets <- read.delim("targets.txt", comment.char = "#")
#Check if targets file is correct
targets
#Create args
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
moduleload(modules(args))

cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="1G") 
# note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
waitForJobs(reg)


