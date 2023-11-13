# load libraries and arguments 
library(rBayesianOptimization)
library(tidyverse)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]
threads = as.numeric(args[4])
n_folds = as.numeric(args[5])

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

ref[1:10, 1:10]
head(labels)

# create n folds 
folds = KFold(labels$label, 
              nfolds = n_folds, 
              stratified = T, 
              seed = 1234)
head(folds)

# Loop through folds and save training and testing data sets 
for (i in 1:n_folds){
  message(paste0('@ SAVING FOLD ', i))
 
  print(head(folds[[i]]))

  # subset test fold
  message('subset test fold')
  test = ref[folds[[i]], ,drop=F]
  test = test %>% rownames_to_column("cell")
  colnames(test)[1] = ""

  # subset true test labels 
  message('subset true test labels')
  test_labels = labels[folds[[i]], ,drop=F]
  test_labels = test_labels %>% rownames_to_column("cell")
  colnames(test_labels)[1] = ""
  
  # subset training data 
  message('subset true test labels')
  train = ref[-folds[[i]], ,drop=F]
  train = train %>% rownames_to_column("cell")
  colnames(train)[1] = "" 
   
  # subset labels for training data
  message('@ subset labels for training data')
  train_labels = labels[-folds[[i]], ,drop=F]
  train_labels = train_labels %>% rownames_to_column("cell")
  colnames(train_labels)[1] = ""

  # save csv files 
  data.table::fwrite(test, paste0(out_path, '/fold', i, '/test.csv'))
  data.table::fwrite(test_labels, paste0(out_path, '/fold', i, '/test_labels.csv'))
  data.table::fwrite(train, paste0(out_path, '/fold', i, '/train.csv'))
  data.table::fwrite(train_labels, paste0(out_path, '/fold', i, '/train_labels.csv'))
}


