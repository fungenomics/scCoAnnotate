# loads needed libraries
library(tidyverse)
library(data.table)
library(Seurat)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

# Script modified from Hussein Lakkis's run_correlation.R script

#--------------- Data -------------------
# read reference matrix and transpose  
message('@ READ REF')
ref = data.table::fread(ref_path, nThread=threads, header=T, data.table=F) %>%
      column_to_rownames('V1') 
message('@ DONE')

# read reference labels
labels = data.table::fread(lab_path, header=T, data.table=F) %>%
         column_to_rownames('V1')


# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(rownames(ref)))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# make Seurat object (transpose ref first)
ref = transposeBigData(ref, blocksize = 10000)
seurat = CreateSeuratObject(counts = ref, meta.data = labels)

# Normalize seurat using default "LogNormalize" method
seurat = NormalizeData(seurat)

#------------- Generate Mean Expression Matrix -------------

# Set labels as identity classes
Idents(seurat) = "label"
# Returns averaged expression values for each identity class
# Function uses slot = "data" (normalized data) as default
seurat = AverageExpression(seurat, return.seurat = TRUE)

# get mean expression matrix
ref_mean_mat = GetAssayData(seurat)

# save mean expression matrix 
message('@ SAVE MODEL')
save(ref_mean_mat, file = out_path)
message('@ DONE')

#----------------------------------------

