# load libraries and arguments 
library(data.table)
library(WGCNA)
library(tidyverse)
library(Seurat)
library(glue)
library(scAnnotate)
set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
query_path = args[3]
pred_path = args[4]
threads = as.numeric(args[5])
threshold = as.numeric(args[6])
# path for other outputs (depends on tools)
out_path = dirname(pred_path)

#--------------- Data -------------------

# read reference matrix and transpose 
message('@ READ REF')
### The matrix of the references is transpose since it needs to be normalize 
### with Seurat that expect a genes x cell.
ref <- data.table::fread(ref_path,
                         data.table=F,
                         header=T,
                         nThread=threads) %>%
  column_to_rownames('V1') %>% 
  transposeBigData()

message('@ DONE')

# read reference labels
labels <- data.table::fread(lab_path,
                            data.table=F,
                            header=T,
                            nThread=threads) %>%
  column_to_rownames('V1') 


# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(colnames(ref)))

# throw error if order is not the same 
if(!order){
  stop("@ Order of cells in reference and labels do not match")
}

# read query matrix and transpose 
message('@ READ QUERY')
### The matrix of the references is transpose since it needs to be normalize 
### with Seurat that expect a genes x cell.
query <- data.table::fread(query_path,
                           data.table=F,
                           header=T,
                           nThread=threads) %>%
  column_to_rownames('V1') %>% 
  WGCNA::transposeBigData() 

message('@ DONE')

### Query and References should have the same gene features
common.genes <- intersect(rownames(query),rownames(ref))
query <- query[common.genes,]
ref   <- ref[common.genes,]
### Prepare the reference
## Normalization
ref <- NormalizeData(ref) %>% as.data.frame() %>% transposeBigData()
## The label should be in the first column
ref <- cbind(labels,ref)

## Prepare the query
## Normalization
query <- NormalizeData(query) %>% as.data.frame() %>% transposeBigData()

#------------- Train + Predict scAnnotate -------------
## Auto to automaticly define the correction method needed.
pred <- scAnnotate(train=ref,
                   test=query,
                   distribution="normal", #Default
                   correction ="auto",
                   screening = "wilcox",
                   threshold=threshold,
                   lognormalized=TRUE)

pred_labels <- data.frame(cell = rownames(query),
                          scAnnotate = pred)

message('@ WRITTING PREDICTIONS')
data.table::fwrite(pred_labels,
                   file = pred_path,
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')
