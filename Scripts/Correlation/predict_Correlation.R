# loads needed libraries
library(tidyverse)
library(data.table)
library(Seurat)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])

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
query = transposeBigData(query, blocksize = 10000)
seurat_query = CreateSeuratObject(counts = query)

# Normalize seurat using default "LogNormalize" method
seurat_query = NormalizeData(seurat_query)

# Get the normalized expression matrix from the query
query_mat = GetAssayData(seurat_query)

# ------------ Load Correlation Function --------------
# Function to label by Spearman correlation generated by Selin Jessa and Marie Coutlier
# Path to original file: /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/stable/code/scripts/predict_celltype_cor.R


label_correlation <- function(test_expr_mat,
                              ref_expr_mat,
                              threshold_common_genes = 0.5) {
  
  rownames(test_expr_mat) <- toupper(rownames(test_expr_mat))
  rownames(ref_expr_mat) <- toupper(rownames(ref_expr_mat))
  
  # Testing how many genes are in common and stopping if not enough
  common_genes <- intersect(rownames(test_expr_mat), rownames(ref_expr_mat))
  prop_common <- length(common_genes) / length(rownames(test_expr_mat))
  message("@@ ", round(prop_common*100, digits = 2), "% of test dataset genes are in the reference dataset")
  
  if (prop_common < threshold_common_genes) stop("Proportion of common genes below threshold.")
  
  # Reducing matrices to common subset
  mat1 <- as.matrix(test_expr_mat[common_genes, ])
  mat2 <- as.matrix(ref_expr_mat[common_genes, ])
  
  # sanity check
  nrow(mat1) == nrow(mat2)
  
  # Computing correlations
  cor_matrix <- cor(mat1, mat2, method = "spearman", use = "complete.obs")
  
  # Getting the best one
  cor_label <- as.data.frame(cor_matrix) %>%
    mutate("cell" = rownames(cor_matrix)) %>%
    gather("celltype", "correlation", -cell) %>%
    group_by(cell) %>%
    top_n(1, correlation) %>%
    dplyr::select(cell, celltype, correlation) %>% 
    arrange(cell)
  
  # Returning the results
  return(cor_label)
  
}

#----------- Predict Correlation ------------

# predict labels 
message('@ PREDICT LABELS')
# Hussein set the threshold_common_genes = 0.3
# This does not affect current prediction as preprocessing is already subsetting common genes between reference and query
# Therefore, 100% of query dataset genes should be in the reference dataset
# However, this is something that could be explored in the future
predicted = as.data.frame(label_correlation(query_mat,ref_mean_mat,0.3))
message('@ DONE')

# Rename columns 
colnames(predicted) <- c("cell", "Correlation","Correlation_score")

# write prediction 
data.table::fwrite(predicted, file = out_path)

#----------------------------------------