# load libraries and arguments 

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])

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

# CODE FOR ANY DATA MANIPULATION NEEDED BEFORE PREDICTION HERE 

#----------- Predict <tool_name> --------

# CODE FOR PREDICTING HERE  
message('@ PREDICT LABELS')


message('@ DONE')

# write prediction 
data.table::fwrite(<pred_table>, file = pred_path)

#----------------------------------------
