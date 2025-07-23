#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# arg[1]: full path to input file (output of getDifferentialPeaksReplicates)
# arg[2]: full path to new output file
# arg[3]: full path to jaspar motif file
# arg[4]: index (numeric) of the column in the input file storing the change in accessibility (log2 FC)

# read inputs
input_file = args[1]
output_file = args[2]
matrix_file = args[3]
logfc_col = as.numeric(args[4])
n_jobs = 16

# imports
library(monaLisa)
library(TFBSTools)
library(readr)
library(genomation)
library(BSgenome.Hsapiens.UCSC.hg38)

# read motifs
pwms = readJASPARMatrix(matrix_file, matrixClass="PWMProb")

# read input data
print("reading input data")
gr = readGeneric(input_file, chr=2, start=3, end=4, strand=5, meta.cols=list(logfc=logfc_col), keep.all.metadata=FALSE, zero.based=FALSE, header=TRUE, skip=0)

# find motif hits in input sequences
print("finding motif hits...")
peakSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
suppressWarnings({
	hits <- findMotifHits(query = pwms, subject = peakSeq, min.score = "85", BPPARAM = BiocParallel::MulticoreParam(n_jobs))
})

# get relevant motifs
TFBSmatrix <- unclass(table(factor(seqnames(hits), levels = seqlevels(hits)), factor(hits$pwmname, levels = name(pwms))))
zero_TF <- colSums(TFBSmatrix) == 0
TFBSmatrix <- TFBSmatrix[, !zero_TF]
print(paste(sum(zero_TF), "motifs were removed"))

# calculate and add G+C and CpG obs/expected to predictor matrix
print("adding GC features to predictor matrix")
fMono <- oligonucleotideFrequency(peakSeq, width = 1L, as.prob = TRUE)
fDi <- oligonucleotideFrequency(peakSeq, width = 2L, as.prob = TRUE)
fracGC <- fMono[, "G"] + fMono[, "C"]
oeCpG <- (fDi[, "CG"] + 0.01) / (fMono[, "G"] * fMono[, "C"] + 0.01)
TFBSmatrix <- cbind(fracGC, oeCpG, TFBSmatrix)

# perform regression
print("performing regression...")
set.seed(123)
se <- randLassoStabSel(x = TFBSmatrix, y = gr$logfc, cutoff = 0.8)

# save selected motif names
print(paste("saving results to", output_file))
selected_motifs = as.vector(colnames(se)[se$selected])
readr::write_lines(selected_motifs, output_file)
print("done")