# the commands input on marcc as an attempt to run MULTI-seq
# first step was to extract the list of cell IDs from the cellranger output
# fastq files I exracted from SRR8890637.sra using fastq-dump from sra toolkit -> 2 fastq files for R1 and R2
# download and extract reference genome from 10x website
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
# unpacking the reference will produce a folder, use this folder name in the cellranger command:
## cellranger count --id=SRR8890637 --transcriptome=refdata-cellranger-GRCh38-3.0.0 --fastqs=fastqs/ --sample=SRR8890637 --localcores=24 --expect-cells=10000
#cell ranger took 6+ hours to complete

# open up an R session on a compute node
## module load R
## module load gcc/5.5.0
# reload R after loading gcc, which allows us to more easily complie and install local R packages on marcc
# install and load all of the following
library(dplyr)
library(ggplot2)
library(Seurat)
library(deMULTIplex)
library(ShortRead)
library(KernSmooth)
library(reshape2)
library(stringdist)
library(Rtsne)

# use seurat to read the cellranger matrix output
cell.mat <- Read10X(data.dir="/home-2/jso4@jhu.edu/scratch/SRR8890637/SRR8890637-cellranger-out/outs/filtered_gene_bc_matrices/GRCh38")
# the column names of the 10x cell ranger matrix are the sequences of the cell IDs
cell.IDs <- data.frame(cell = colnames(cell.mat))
# record sample cell IDs in another file
write.table("/home-2/jso4@jhu.edu/scratch/SRR8890637/SRR8890637_cellIDs.txt", sep = "\t", col.names=F, row.names=F, quote=F)
cell.id.vec <- as.character(cell.IDs$cell)

# import the list of MULTI-seq barcodes which were provided in an excel spreadsheet
# I uploaded them to marcc as a .txt file
bar.ref <- read.delim("multiBC.txt", sep="\n", header=F)
# make it into a character vector
bar.ref <- as.character(bar.ref$V1)

# run MULTI-seq preProcessing function
# specify UMI as c(17,26) for V2 chemistry, c(17,28) for V3
readTable <- MULTIseq.preProcess(R1="SRR8890637/fastqs/SRR8890637_S1_L001_R1_001.fastq.gz", R2="SRR8890637/fastqs/SRR8890637_S1_L001_R2_001.fastq.gz", cellIDs=cell.id.vec, cell=c(1,16), umi=c(17,26), tag=c(1,8))
# save output as a .csv
write.csv(readTable, file="SRR8890637_readTable.csv", quote=F)

# perform MULTI-seq sample barcode alignment
bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
# save result
write.csv(bar.table, file="SRR8890637_barcode_counts.csv", quote=F)

# Visualize barcode space
bar.tsne <- barTSNE(bar.table[,1:96])

# paranoia, save all the things
save.tsne <- bar.table
save.tsne$TSNE1 <- bar.tsne$TSNE1
save.tsne$TSNE2 <- bar.tsne$TSNE2
write.csv(save.tsne, file="SRR8890637_barcode_tsne.csv", quote=F)

# since I can't get the for loops to work on marcc,
## I'm going to copy the output .csv files to my desktop and try to do the analysis in R studio
