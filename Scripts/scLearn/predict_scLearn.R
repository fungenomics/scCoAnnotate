# load libraries and arguments 
library(tidyverse)
library(M3Drop)
library(scLearn)
library(glue)
set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
# species = args[5] # species="Hs" for homo sapiens or species="Mm" for mus musculus.

# path for other outputs (depends on tools)
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

# transpose query data so so cells are on columns
# Normalizing with Seurat, in order to not filter any cell.
query <- query %>% WGCNA::transposeBigData() %>% Seurat::NormalizeData()

#----------- Predict scLearn --------

# CODE FOR PREDICTING HERE  
message('@ PREDICT LABELS')

# From documentation: "Assignment with trained model above. To get a less strict result for "unassigned" cells, you can decrease "diff" and "vote_rate". If you are sure that the cell type of query cells must be in the reference dataset, you can set "threshold_use" as FALSE. It means you don't want to use the thresholds learned by scLearn."
scLearn_predict_result <- scLearn_cell_assignment(scLearn_model_learning_result = scLearn,
                                                  expression_profile_query = query,
                                                  diff = 0.05,
                                                  threshold_use = TRUE,
                                                  vote_rate = 0.6)

message('@ DONE')


pred_labs = data.frame(cell = scLearn_predict_result$Query_cell_id,
                       scLearn = scLearn_predict_result$Predict_cell_type)

# write prediction 
message('@ WRITE PREDICTIONS')
data.table::fwrite(pred_labs, 
                   file = pred_path,
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')
#----------------------------------------
### Add the entire output with the Scores in Additional_information
message('@ WRITE TABLE OUTPUT WITH ALL PROBABILITIES')
data.table::fwrite(scLearn_predict_result, 
                   file = glue('{out_path}/table_with_probabilities.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')