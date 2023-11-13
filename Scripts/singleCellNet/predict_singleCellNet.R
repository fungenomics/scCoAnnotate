library(tidyverse)
library(Seurat)
library(singleCellNet)
library(tidyverse)
library(data.table)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

#--------------- Data -------------------

# read query matrix and transpose 
message('@ READ QUERY')
query = data.table::fread(sample_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# Get cell names
cellnames = row.names(query)

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# Transpose query 
query = transposeBigData(query, blocksize = 10000)

#----------- Predict singleCellNet ------------

# predict labels 
message('@ PREDICT LABELS')

# Default nrand = 50, Hussein had used nrand = 0
pred = scn_predict(class_info[['cnProc']], query, nrand = 0)
message('@ DONE')

# classify cells 
query_res = assign_cate(classRes = pred, sampTab = data.frame(row.names = cellnames), cThresh = 0.5) 

pred_labs = data.frame(cell = rownames(query_res),
	               singleCellNet = query_res$category)

# write prediction 
data.table::fwrite(pred_labs, file = out_path)

#----------------------------------------
