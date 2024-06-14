# load libraries and arguments 
library(tidyverse)
library(data.table)
library(Seurat)
library(WGCNA)

set.seed(1234)

#---------- Parameters -------------------
args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])
nPC_computed = as.numeric(args[5])

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
seurat_ref = CreateSeuratObject(counts = ref,
                                meta.data = labels)


# if('batch' %in% colnames(seurat_ref@meta.data)){
#   seurat_ref[["RNA"]] <- base::split(seurat_ref[['RNA']],
#                                f = seurat_ref$batch)
# }
# Normalize seurat using default "LogNormalize" method
seurat_ref = NormalizeData(seurat_ref)
seurat_ref = FindVariableFeatures(seurat_ref)
seurat_ref = ScaleData(seurat_ref)
seurat_ref = RunPCA(seurat_ref,
                    npcs = nPC_computed)

# ## If the batch on the reference were specified, performs an integration
# seurat_ref <- IntegrateLayers(object = seurat_ref,
#                               method = FastMNNIntegration,
#                               orig.reduction = "pca",
#                               new.reduction = "integrated.cca",
#                               verbose = TRUE)


# save the pre-processed object
message('@ SAVE MODEL')
save(seurat_ref,
     file = out_path)
message('@ DONE')

#----------------------------------------

