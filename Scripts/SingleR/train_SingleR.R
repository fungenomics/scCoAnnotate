# load libraries
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

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
order = all(rownames(labels..) == rownames(ref..))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# make SingleCellExperiment object
ref = SingleCellExperiment(assays = list(counts = t(ref)))

# log normalize reference 
message('@ NORMALIZE REF')
ref = scuttle::logNormCounts(ref)
message('@ DONE')

#------------- Train SingleR -------------

# train SingleR
message('@ TRAINING MODEL')
singler = trainSingleR(ref, labels=labels$label, num.threads = threads, assay.type = "logcounts")
message('@ DONE')

# save trained model 
message('@ SAVE MODEL')
save(singler, file = paste0(out_path, '/model_SingleR.Rda'))
message('@ DONE')

#----------------------------------------
