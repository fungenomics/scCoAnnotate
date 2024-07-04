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

pred_matrix <- Scibet_model(query,
                            result = 'table')
rownames(pred_matrix) <- rownames(query)

##Get the max value score in the probabiliy matrix for each cell is the same as run the Scibet_model with result = 'list'
## I take the max value for each row
pred <- apply(pred_matrix,
              1,
              function(x){
                names(which.max(x))
                }
              )
#Since is a named vector I use the names (cells) and the values (labels)
pred_labels <- data.frame(cell = names(pred),
                          SciBet = as.character(pred))

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
pred_matrix <- pred_matrix %>% as.data.frame() %>% tibble::rownames_to_column(" ")

message('@ SAVE PRED MATRIX')
data.table::fwrite(pred_matrix,
                   file = glue('{out_path}/SciBet_pred_score.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads
)
message('@ DONE')
#----------------------------------------
