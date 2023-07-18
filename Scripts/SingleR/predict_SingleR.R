# load libraries
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

#--------------- Data -------------------

# read query matrix and transpose 
message('@ READ QUERY')
query = data.table::fread(sample_path, nThread=threads, header=T, data.table=F, nrow = 500) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# Make SingleCellExperiment object from query
query = SingleCellExperiment(assays = list(counts = query))

# Log normalize query counts 
message('@ NORMALIZE QUERY')
query = scuttle::logNormCounts(query)
message('@ DONE')
print(query)


#----------- Predict SingleR ------------

# predict labels 
message('@ PREDICT LABELS')
pred = classifySingleR(query, singler, assay.type = "logcounts")
message('@ DONE')

save(pred, file = paste0(out_path, '/SingleR_pred.Rda'))

pred_labs = data.frame(SingleR = pred$labels,
	               row.names = rownames(pred))

print(pred_labs)

# write prediction 
data.table::fwrite(pred_labs, file = paste0(out_path, '/SingleR_pred.csv'))

#----------------------------------------

