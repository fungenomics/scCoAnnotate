# scClassify prediction script
# Bhavyaa Chandarana, July 2023

# load libraries and arguments 
library(data.table) 
library(Seurat)
library(WGCNA)
library(scClassify)
library(dplyr)
library(tibble)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])

# path for other outputs (depends on tools)
out_path = dirname(model_path)

#--------------- Data -------------------

# read query matrix 
message('@ READ QUERY')
query <- data.table::fread(query_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
# object name "scClassify"
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# transpose (put cells in columns) for Seurat normalization and scClassify, normalize 
query <- query %>% WGCNA::transposeBigData() %>% Seurat::NormalizeData()

#----------- Predict scClassify --------

# specify parallelization configuration depending on number of threads
if(threads > 1){

    bpparam <- BiocParallel::MulticoreParam(workers = threads)
    parallel <- TRUE

} else {
   
    bpparam <- BiocParallel::SerialParam()
    parallel <- FALSE

}

# run prediction
message('@ PREDICTING QUERY')
pred <- predict_scClassify(
    exprsMat_test = as.matrix(query),
    trainRes = scClassify,
    parallel = parallel,
    BPPARAM = bpparam
)

# extract predictions table and format
message('@ FORMATTING PREDICTIONS')
pred_table <- pred$pearson_WKNN_limma$predRes %>% as.data.frame() %>% rownames_to_column()
colnames(pred_table)[1] <- "cell"
colnames(pred_table)[2] <- "scClassify"

# write prediction 
message('@ SAVE PREDICTIONS')
data.table::fwrite(pred_table, file = pred_path,
                    row.names = F,
                    col.names = T,
                    sep = ",",
                    nThread = threads)
message('@ DONE')

#----------------------------------------
