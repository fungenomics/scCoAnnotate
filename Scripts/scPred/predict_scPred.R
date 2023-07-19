# load libraries and arguments 

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

#--------------- Data -------------------

# read query matrix 
message('@ READ QUERY')
query = data.table::fread(sample_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# transpose and create seurat object 
query = t(query)
query = CreateSeuratObject(query, row.names = colnames(query))

# normalize query 
query = query %>% 
  NormalizeData() 

#----------- Predict scPred -------------

# predict cells 
query = scPredict(query, scpred)

pred_labs = data.frame(cell = rownames(query),
                       scPred = query$scpred)

# write prediction 
data.table::fwrite(pred_labs, file = out_path)

#----------------------------------------