# load libraries and arguments 
library(tidyverse)
library(data.table)
library(symphony)
library(WGCNA)
library(Seurat)
library(glue)
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

ref = transposeBigData(ref, blocksize = 10000) %>% Seurat::as.sparse()

# https://github.com/immunogenomics/symphony/issues/9
batch <- NULL
if('batch' %in% colnames(labels)){
  batch <- 'batch'
  labels$batch <- paste0('ref_',labels$batch) 
  message('@Running with batches')
}


## I use all the default values provided here: https://github.com/immunogenomics/symphony?tab=readme-ov-file
## If vars is null no batch correction is run
## The default nPC_computed is 50 to match with the seurat pieline
  
# Build reference
symphony_ref = symphony::buildReference(
  ref,                       # reference expression (genes by cells)
  labels,                    # reference metadata (cells x attributes)
  vars = batch,              # variable(s) to integrate over
  K = 100,                   # number of Harmony soft clusters
  verbose = TRUE,            # display verbose output
  do_umap = FALSE,           # run UMAP and save UMAP model to file
  do_normalize = TRUE,       # perform log(CP10k) normalization on reference expression
  vargenes_method = 'vst',   # variable gene selection method: 'vst' or 'mvp'
  vargenes_groups = batch,   # metadata column specifying groups for variable gene selection within each group
  topn = 2000,               # number of variable genes (per group)
  theta = 2,                 # Harmony parameter(s) for diversity term
  d = nPC_computed,          # number of dimensions for PCA
  additional_genes = NULL    # vector of any additional genes to force include
)

# save the pre-processed object
message('@ SAVE MODEL')
save(symphony_ref,
     file = out_path)
message('@ DONE')

#----------------------------------------
