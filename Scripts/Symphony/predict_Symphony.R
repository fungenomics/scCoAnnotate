# load libraries and arguments 
library(tidyverse)
library(data.table)
library(Seurat)
library(WGCNA)
library(symphony)

set.seed(1234)

#---------- Parameters -------------------
args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
# nPC_used = as.numeric(args[5])

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

# Transpose query 
query = transposeBigData(query, blocksize = 10000) %>% Seurat::as.sparse()

# Map query
# It could be fixed the batch in the query, but this is not implemented yet
# for the metadata_query I create an empty data.frame because it's needed afer for the kNN prediction
query = mapQuery(exp_query = query,     # query gene expression (genes x cells) 
                 metadata_query = data.frame(row.names = colnames(query)), # query metadata (cells x attributes)
                 ref_obj =  symphony_ref,  # Symphony reference object
                 vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                 do_normalize = TRUE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                 do_umap = FALSE)        # project query cells into reference UMAP
# extract predictions table and format

query = knnPredict(query_obj = query,
                   ref_obj = symphony_ref,
                   train_labels = symphony_ref$meta_data$label,
                   k = 5,
                   # save_as = 'predicted_labels',
                   confidence = F
                   )

message('@ FORMATTING PREDICTIONS')
pred_labs <- data.frame(cell = rownames(query$meta_data),
                         Symphony = query$meta_data$cell_type_pred_knn
                         )
# write prediction 
message('@ SAVE PREDICTIONS')
data.table::fwrite(pred_labs,
                   file = pred_path,
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')

#----------------------------------------

# output binary matrix
# message('@ WRITE TABLE WITH BINARY OUTPUT')
# pred_labs$Symphony <- factor(pred_labs$Symphony,
#                              levels = unique(symphony_ref$meta_data$label)) #set all the posible values

pred_labs = pred_labs %>% 
  mutate(prob = 1) %>% 
  pivot_wider(names_from = Symphony, 
              values_from = prob, 
              values_fill = 0) %>% as.data.frame
              # names_expand = TRUE) %>% as.data.frame() 
colnames(pred_labs)[1] <- " "

data.table::fwrite(pred_labs, 
                   file = paste0(out_path, '/Symphony_pred_score.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')

## Save the entire query output
save(query,
     file = paste0(out_path, '/query_results.Rds')
)