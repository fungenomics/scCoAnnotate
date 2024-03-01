# scClassify prediction script
# Bhavyaa Chandarana, July 2023

# load libraries and arguments 
library(data.table) 
library(Seurat)
library(scClassify)
library(tidyverse)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
algorithm = args[5]
similarity = args[6]
prob_threshold = as.numeric(args[7])
cor_threshold_static = as.numeric(args[8])
cor_threshold_high = as.numeric(args[9])

# path for other outputs (depends on tools)
out_path = dirname(pred_path)

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
query <- query %>% transposeBigData() %>% Seurat::NormalizeData()

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
    algorithm = 
    trainRes = scClassify,
    parallel = parallel,
    BPPARAM = bpparam
)

print(head(pred))

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

# output binary matrix
message('@ WRITE TABLE WITH BINARY OUTPUT')
pred_table = pred_table %>% 
            mutate(prob = 1) %>% 
            pivot_wider(names_from = scClassify, 
                        values_from = prob, 
                        values_fill = 0)

names(pred_table)[1] = ""

data.table::fwrite(pred_table, 
                   file = paste0(out_path, '/scClassify_pred_score.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')
