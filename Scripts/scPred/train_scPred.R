# load libraries and arguments 
library(scPred)
library(Seurat)
library(tidyverse)
library(doParallel)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
model_path = args[3]
threads = as.numeric(args[4])
model_type = args[5]

# path for other outputs (depends on tools)
out_path = dirname(model_path)

#--------------- Data -------------------

# read reference matrix 
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

# transpose and create seurat object with the labels as meta data 
ref = transposeBigData(ref)
ref = CreateSeuratObject(ref, meta.data = labels)

# normalize reference
ref = ref %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

# create scPred object stored in the @misc slot of reference 
ref = getFeatureSpace(ref, "label")

#------------- Train scPred -------------

# train model (parallelized) 
cl = makePSOCKcluster(threads)
registerDoParallel(cl)
scpred = trainModel(ref, model = model_type, allowParallel = T)
stopCluster(cl)

# print model info 
get_scpred(scpred)

# save trained model
message('@ SAVE MODEL')
save(scpred, file = model_path)
message('@ DONE')

# Plot prob 
pdf(paste0(out_path, '/qc_plots.pdf'), width=10, height=10)
plot_probabilities(scpred)
dev.off()

#---------------------------------------------
