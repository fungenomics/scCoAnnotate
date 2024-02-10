# load libraries
library(tidyverse)
library(data.table)
library(Seurat)
library(singleCellNet)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])
nTrees = as.numeric(args[5])

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

#------------- Train singleCellNet -------------

# Split reference data for training and assessment
# Default ncells = 50, Hussein had used ncells = 80 in scCoAnnotate V1
stList = splitCommon(sampTab = seurat@meta.data, ncells = 50, dLevel = "label")
# Get the downsampled list
stTrain = stList[[1]]
# Get corresponding 
expTrain = as.matrix(GetAssayData(seurat))[,row.names(stTrain)]

# Train singleCellNet
# Default uses nTopGenes = 10, nTrees = 1000
# Hussein had used nTopGenes = 12, nTrees = 350 in scCoAnnotate V1 and nTrees = 500 in his thesis
message('@ TRAINING MODEL')
class_info = scn_train(stTrain = stTrain, 
                        expTrain = expTrain, 
                        nTopGenes = 10, 
                        nRand = 70, 
                        nTrees = nTrees, 
                        nTopGenePairs = 25, 
                        dLevel = "label")
message('@ DONE')

# save trained model 
message('@ SAVE MODEL')
save(class_info, file = out_path)
message('@ DONE')

#----------------------------------------




