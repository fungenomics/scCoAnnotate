# load libraries and arguments 
library(data.table)
library(scibet)
library(WGCNA)
library(tidyverse)
library(Seurat)
library(glue)

set.seed(1234)

#---------- Parameters -------------------
args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
# path for other outputs (depends on tools)
out_path = dirname(pred_path)

#--------------- Data -------------------

# read query matrix and transpose 
message('@ READ QUERY')

# The matrix of the references is transpose since it needs to be normalize 
# with Seurat that expect a genes x cell.
query <- data.table::fread(query_path,
                          data.table=F,
                          header=T,
                          nThread=threads) %>%
  column_to_rownames('V1') %>% 
  WGCNA::transposeBigData() 

message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# The matrix here is transposed since SciBet expect a cell x gene matrix.
query <- Seurat::NormalizeData(query) %>% as.data.frame() %>% WGCNA::transposeBigData()

#----------- Predict SciBet --------

pred <- Scibet_model(query,
                     result = 'list')

pred_labels <- data.frame(cell = rownames(query),
                          SciBet = pred)

message('@ WRITTING PREDICTIONS')
data.table::fwrite(pred_labels,
                   file = pred_path,
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')

#------------- Other outputs --------------

# I tested and the results are the same when you run the same model twice,
# so I run it again to obtain the prob matrix.
pred_matrix <- Scibet_model(query,
                            result = 'table')
rownames(pred_matrix) <- rownames(query)
pred_matrix <- pred_matrix %>% as.data.frame() %>% tibble::rownames_to_column(" ")

message('@ SAVE PRED MATRIX')
data.table::fwrite(pred_matrix,
                   file = glue('{out_path}/prob_matrix.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads
)
message('@ DONE')
#----------------------------------------
