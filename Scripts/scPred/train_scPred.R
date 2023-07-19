# load libraries and arguments 
library(scPred)
library(Seurat)
library(tidyverse)
library(doParallel)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])
model_dir = args[5]
model = args[6]

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
ref = t(ref)
ref = CreateSeuratObject(ref, row.names = colnames(ref), meta.data = labels)

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
scpred = trainModel(ref, model = model, allowParallel = T)
stopCluster(cl)

# print model info 
get_scpred(scpred)

# Plot prob 
pdf(paste0(model_dir, 'qc_plots.pdf'), width=10, height=10)
plot_probabilities(scpred)
dev.off()

# save trained model 
message('@ SAVE MODEL')
save(scpred, file = out_path)
message('@ DONE')

#---------------------------------------------
