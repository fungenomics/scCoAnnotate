# load libraries and arguments 
library(data.table)
library(WGCNA)
library(tidyverse)
library(glue)
library(Seurat)
library(scID)
library(MAST)
set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
query_path = args[3]
pred_path = args[4]
threads = as.numeric(args[5])
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

### Prepare the reference
## Normalization
ref <- scID:::counts_to_cpm(counts_gem = ref)

## Labels as named vector
label <- setNames(object = labels$label,nm = rownames(labels))
## Prepare the query
## Normalization
query <- scID:::counts_to_cpm(counts_gem = query)


#----------- Train + Predict scID --------
pred <- scid_multiclass(target_gem = query,
                               reference_gem = ref,
                               reference_clusters = label, 
                               logFC = 0.5, #Default
                               only_pos = FALSE, #Default
                               normalize_reference = FALSE, #I already normalized the reference
                               estimate_weights_from_target = FALSE #Default
                               )

pred_labels <- data.frame(cell = names(pred$labels),
                          scID = as.character(pred$labels)
                          )

message('@ WRITTING PREDICTIONS')
data.table::fwrite(pred_labels,
                   file = pred_path,
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')

#------------- Other outputs --------------
#Save the entire output
message('@ SAVE scID OUTPUT')
save(pred,
     file =  glue('{out_path}/scID_output.Rdata')
     )
message('@ DONE')