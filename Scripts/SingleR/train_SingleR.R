# load libraries
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
model_path = args[3]
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
order = all(as.character(rownames(labels)) == as.character(rownames(ref)))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# make SingleCellExperiment object (transpose ref first)
ref = transposeBigData(ref)
ref = SingleCellExperiment(assays = list(counts = ref))

# log normalize reference 
message('@ NORMALIZE REF')
ref = scuttle::logNormCounts(ref)
message('@ DONE')

#------------- Train SingleR -------------

# train SingleR
message('@ TRAINING MODEL')
singler = trainSingleR(ref, 
                       labels=labels$label, 
		       assay.type = "logcounts",
		       de.method = "wilcox",
		       de.n = 10)
message('@ DONE')

# save trained model 
message('@ SAVE MODEL')
save(singler, file = model_path)
message('@ DONE')

#----------------------------------------
