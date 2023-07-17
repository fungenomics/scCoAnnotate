# load libraries
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

#--------------- Data --------------

# read reference matrix and transpose 
message('@ READ REF')
ref = data.table::fread(ref_path, data.table=F, header=T, nThread=threads) %>%
      as.data.frame() %>% 
      column_to_rownames('V1') %>%
      t()
message('@ DONE')

# Make SingleCellExperiment object
ref = SingleCellExperiment(assays = list(counts = ref))

# Log normalize reference 
message('@ NORMALIZE REF')
ref = scuttle::logNormCounts(ref)
message('@ DONE')

# Read Reference labels
labels = data.table::fread(lab_path)

# check if cell names are in the same order in labels and ref
order = all(labels$cell == colnames(ref))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

#------------- SingleR -------------

# Train SingleR
message('@ TRAINING MODEL')
singler = trainSingleR(ref, labels=labels$label, num.threads = threads)
message('@ DONE')

# save trained model 
message('@ SAVE MODEL')
save(singler, file = paste0(out_path, '/model_SingleR.Rda'))
message('@ DONE')

#-----------------------------------
