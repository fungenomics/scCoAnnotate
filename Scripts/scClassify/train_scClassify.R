# scClassify training script
# Bhavyaa Chandarana, July 2023

# load libraries and arguments
library(data.table) 
library(Seurat)
library(scClassify)
library(tidyverse)
library(ggplot2)
library(glue)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
model_path = args[3]
threads = as.numeric(args[4])

# path for other outputs (depends on tools)
out_path = dirname(model_path)

#--------------- Data -------------------

# read reference matrix 
message('@ READ REF')
ref <- data.table::fread(ref_path, nThread=threads, header=T, data.table=F) %>% 
       column_to_rownames("V1")
message('@ DONE')

# read reference labels
labels <- data.table::fread(lab_path, header=T, data.table=F) %>%
         column_to_rownames('V1')
message('@ DONE')

# check if cell names are in the same order in labels and ref
order <- all(as.character(rownames(labels)) == as.character(rownames(ref)))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# transpose (put cells in columns) for Seurat normalization and scClassify, normalize 
ref <- ref %>% transposeBigData() %>% Seurat::NormalizeData()

#------------- Train scClassify -------------

# specify parallelization configuration depending on number of threads
if(threads > 1){

    bpparam <- BiocParallel::MulticoreParam(workers = threads)
    parallel <- TRUE

} else {
   
    bpparam <- BiocParallel::SerialParam()
    parallel <- FALSE

}

# train scClassify  
message('@ TRAINING MODEL')
scClassify <- train_scClassify(
    exprsMat_train = as.matrix(ref),
    cellTypes_train = labels$label,
    parallel = parallel,
    BPPARAM = bpparam,
    returnList = FALSE # so that cell type tree can be extracted
)

# save trained model 
message('@ SAVE MODEL')
save(scClassify, file = model_path)
message('@ DONE')

#------------- Other outputs ----------------

# save cell type tree produced in training
message('@ SAVE TREE')
tree <- cellTypeTree(scClassify) 
save(tree, file = glue("{out_path}/scClassify_tree.Rda"))
message('@ DONE')

# save plot of same cell type tree
message('@ SAVE TREE PLOT')
tree_plot <- plotCellTypeTree(tree)
ggsave(tree_plot, file = glue("{out_path}/scClassify_tree.png"))
message('@ DONE')

#---------------------------------------------
