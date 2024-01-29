# load libraries and arguments 
library(scPred)
library(Seurat)
library(tidyverse)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])

# get path for other output
out_path = dirname(pred_path)

#--------------- Data -------------------

# read query matrix 
message('@ READ QUERY')
query = data.table::fread(query_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# transpose and create seurat object 
query = transposeBigData(query)
query = CreateSeuratObject(query)

# normalize query 
query = query %>% 
  NormalizeData() 

#----------- Predict scPred -------------

# predict cells 
message('@ PREDICT LABELS')
query = scPredict(query, scpred)
message('@ DONE')

head(colnames(query))
head(query$scpred_prediction)

# scPred chnages - to _minus --> chnage back before saving 
query$scpred_prediction = gsub("_minus", "-", query$scpred_prediction)
pred_labs = data.frame(cell = colnames(query),
                       scPred = query$scpred_prediction)

# write prediction 
data.table::fwrite(pred_labs, file = pred_path)

# save probbability matrix 
prob_mat = query@meta.data %>% rownames_to_column('cell') %>% select(cell | starts_with('scpred'))
colnames(prob_mat) = str_remove(colnames(prob_mat), 'scpred_')
colnames(prob_mat) = gsub("_minus", "-", colnames(prob_mat))
colnames(prob_mat)[1] = ""

# write probability matrix 
data.table::fwrite(prob_mat, file = paste0(out_path, '/scPred_pred_score.csv'))

#----------------------------------------
