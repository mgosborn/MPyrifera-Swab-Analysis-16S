#This was built using https://benjjneb.github.io/dada2/tutorial.html

setwd('/staging/sn1/mgosborn/fulltest/PairedReads/dada2')
print("Set working directory to /staging/sn1/mgosborn/fulltest/PairedReads/dada2")
print("Loading dada2")
library(dada2); packageVersion("dada2")
print("dada2 loaded")

path <- "/staging/sn1/mgosborn/fulltest/PairedReads/dada2" # CHANGE ME to the directory containing the fastq files after unzipping.
print("Files in directory:")
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
print("Extracted sample names")

#print("Plotting quality of fwd reads.")
#pdf("FwdQuality")
#plotQualityProfile(fnFs[])
#dev.off()

#print("Plotting quality of rev reads.")
#pdf("RevQuality")
#plotQualityProfile(fnRs[])
#dev.off()

# Place filtered files in filtered/ subdirectory
print("Placing filtered reads in filtered/ subdirectory")
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(122,130),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On cluster, set multithread = FALSE
print("Reads trimmed and filtered.")
print("Here's a preview:")
head(out)

print("Learning error rates...")
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
print("Errors learned. Now plotting...")

pdf("ErrorRates")
plotErrors(errF, nominalQ=TRUE)
dev.off()

print("Applying core inference algorithm...")
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
#By default, the dada function processes each sample independently. 
#However, pooling information across samples can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples. 
#The dada2 package offers two types of pooling. dada(..., pool=TRUE) performs standard pooled processing, in which all samples are pooled together for sample inference. 
#dada(..., pool="pseudo") performs pseudo-pooling, in which samples are processed independently after sharing information between samples, approximating pooled sample inference in linear time.

print("Merging reads...")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)
# Inspect the merger data.frame from the first sample
print("Reads merged. Here's a preview:")
head(mergers[[1]])

print("Constructing a sequence table...")
seqtab <- makeSequenceTable(mergers)
print("The dimensions are:")
dim(seqtab)
print("Distribution of sequence lengths is:")
table(nchar(getSequences(seqtab)))

print("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
print("The dimensions of new sequence table is:")
dim(seqtab.nochim)
print("Percent of reads remaining:")
sum(seqtab.nochim)/sum(seqtab)
print("Here's a preview of the number of reads that made it through each step in the pipeline:")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

print("Assigning taxonomy...")
taxa <- assignTaxonomy(seqtab.nochim, "/staging/sn1/mgosborn/fulltest/PairedReads/dada2/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
print("Done! Here's a preview:")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

print("No mock community present. Will not evaluate accuracy.")

print("Boom shakalaka!")