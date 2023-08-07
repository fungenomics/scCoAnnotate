# load libraries and arguments 
library(tidyverse)
library(M3Drop)
library(scLearn)
library(Seurat)
library(glue)
set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
model_path = args[3]
threads = as.numeric(args[4])
#species = args[5] # species="Hs" for homo sapiens or species="Mm" for mus musculus.

# path for other outputs (depends on tools)
out_path = dirname(model_path)

#--------------- Data -------------------

# read reference matrix 
message('@ READ REF')
ref = data.table::fread(ref_path, nThread=threads, header=T, data.table=F) %>%
      column_to_rownames('V1')
message('@ DONE')

# read reference labels
labels = data.table::fread(lab_path, header=T, data.table=F) %>%
         column_to_rownames('V1')

# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(rownames(ref)))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# change labels format to named list (names = cell IDs, values = labels)
labels <- setNames(nm=rownames(labels),
                   object = labels$label)



# Preprocessing 
# transpose reference so cells are on columns
## Because we don't want to filter any cell, we will do only the normalization as they do it.
### They do manually the Seurat normalization.
# ref <-apply(ref,2,function(x){x/(sum(x)/10000)})
# ref <- log(ref+1)

ref <- ref %>% WGCNA::transposeBigData() %>% Seurat::NormalizeData()

### Select genes
high_varGene_names <- Feature_selection_M3Drop(expression_profile = ref,
                                               log_normalized = T #True because it was previously normalized
                                               )
#------------- Train scLearn -------------

# From documentation: "training the model. To improve the accuracy for "unassigned" cell, you can increase "bootstrap_times", but it will takes longer time. The default value of "bootstrap_times" is 10."
scLearn <- scLearn_model_learning(high_varGene_names = high_varGene_names,
                                  expression_profile = as.matrix(ref), #It need a matrix not a sparse matrix.
                                  sample_information_cellType = labels,
                                  bootstrap_times = 10 #Default
                                  ) 

# save trained model 
message('@ SAVE MODEL')
save(scLearn,
     file = model_path)
message('@ DONE')

#---------------------------------------------

data.table::fwrite(data.frame(Selected_genes = scLearn$high_varGene_names), 
                   file = glue('{out_path}/selected_genes.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads)
