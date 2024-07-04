# load libraries
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])

# get path for other output
out_path = dirname(pred_path)

#--------------- Data -------------------

# read query matrix and transpose 
message('@ READ QUERY')
query = data.table::fread(sample_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# Make SingleCellExperiment object from query (transpose query first)
query = transposeBigData(query)
query = SingleCellExperiment(assays = list(counts = query))

# Log normalize query counts 
message('@ NORMALIZE QUERY')
query = scuttle::logNormCounts(query)
message('@ DONE')

#----------- Predict SingleR ------------

# predict labels 
message('@ PREDICT LABELS')
pred = classifySingleR(query, singler, assay.type = "logcounts")
message('@ DONE')

head(pred)

pred_labs = data.frame(cell = rownames(pred),
	               SingleR = pred$labels)

# write prediction 
data.table::fwrite(pred_labs, file = pred_path)

# save correlaion matrix 
cor = pred$scores %>% as.data.frame() %>% `rownames<-`(rownames(pred)) %>% rownames_to_column('cell')
colnames(cor)[1] = ""

data.table::fwrite(cor, file = paste0(out_path, '/SingleR_pred_score.csv'))

#----------------------------------------

