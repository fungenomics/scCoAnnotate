# load libraries and arguments 
library(tibble)
library(data.table)
library(Seurat)
library(WGCNA)

set.seed(1234)

#---------- Parameters -------------------
args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
nPC_used = as.numeric(args[5])
integration_method <- args[6]
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

pca.use <- 'pca'
if('integrated.dim' %in% names(seurat_ref@reductions)){
  message('@Using integration')
  pca.use <- 'integrated.dim'
}

# Transpose query 
query = transposeBigData(query, blocksize = 10000)
seurat_query = CreateSeuratObject(counts = query)

# Normalize seurat using default "LogNormalize" method
seurat_query = NormalizeData(seurat_query)
seurat_query.anchors <- FindTransferAnchors(reference = seurat_ref,
                                            query = seurat_query,
                                            dims = 1:nPC_used,
                                            reference.reduction = pca.use)

predictions <- TransferData(anchorset = seurat_query.anchors,
                            refdata = seurat_ref$label,
                            dims = 1:nPC_used)

# extract predictions table and format
message('@ FORMATTING PREDICTIONS')
pred_table <- data.frame(cell = rownames(predictions),
                         Seurat = predictions$predicted.id)

# write prediction 
message('@ SAVE PREDICTIONS')
data.table::fwrite(pred_table,
                   file = pred_path,
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
message('@ DONE')

# write score matrix 
prob_mat <- predictions
colnames(prob_mat) <- gsub(x = colnames(predictions),
                              pattern = "^prediction.score.",
                              replacement = "")
prob_mat <- prob_mat[,-which(colnames(prob_mat) %in% c("max","predicted.id"))]
prob_mat <- prob_mat %>% tibble::rownames_to_column(" ")
## The sum per cell is 1, so we can consider them as a probability matrix
data.table::fwrite(prob_mat, 
                   file = paste0(out_path, '/Seurat_pred_score.csv'))