# load libraries and arguments 
library(rBayesianOptimization)
library(tidyverse)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
out_path = args[3]

threads = as.numeric(args[4])
if(is.na(threads)){
  stop("The number threads specified is not a numeric value")
}

n_folds = as.numeric(args[5])
if(is.na(n_folds)){
  stop("The number of folds specified is not a numeric value")
}

min_cells = as.numeric(args[6])
if(is.na(min_cells)){
  stop("The minimum number of cells specified is not a numeric value")
}

downsample_value = as.numeric(args[7]) 
if(is.na(downsample_value)){
  stop("The downsample value specified is not a numeric value")
}
print(class(downsample_value))
downsample_stratified = as.logical(args[8])
if(is.na(downsample_stratified)){
  stop("The downsample stratified specified is not a logical value")
}
downsample_stratified = if(downsample_stratified) "label" else NULL

print(downsample_stratified)
print(class(downsample_stratified))
#--------------- Data -------------------

# read reference matrix 
message('@ READ REF')
ref = data.table::fread(ref_path, nThread=threads, header=T, data.table=F) %>%
      column_to_rownames('V1')
message('@ DONE')

# read reference labels
labels = data.table::fread(lab_path, header=T, data.table=F) 

if(downsample_value != 0){
  if(downsample_value >= 1){
    labels = labels %>%  group_by(across(all_of(downsample_stratified))) %>% 
      dplyr::slice_sample(n = downsample_value,replace = F)
  } else{
    labels = labels %>% group_by(across(all_of(downsample_stratified))) %>% 
      dplyr::slice_sample(prop = downsample_value,replace = F) 
  }
  ref = ref[labels$V1,]
}
#do the convertion to rownames after since the column is needed to downsampling
labels = labels %>% column_to_rownames('V1')


# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(rownames(ref)))

if(min_cells > 0){
  rmv_labels = names(which(table(labels$label) < min_cells))
  labels = labels %>% filter((!label %in% rmv_labels))
  message(paste0(paste0(rmv_labels,collapse = '-'),' classes were remove because of lower number of cells (< ',as.character(min_cells),')'))
  #filtering the cells from the filtered classes
  ref = ref[rownames(labels),]
}

data.table::fwrite(data.frame(cells= rownames(labels),
                              label= labels$label),
                   file = paste0(out_path,'/downsampled_reference_labels.csv'), sep = ',')

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# ref[1:10, 1:10]
# head(labels)

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

lab = data.frame(label = unique(labels$label))
data.table::fwrite(lab,
                   file = paste0(out_path, '/ontology/ontology.csv'),
                   sep = ',')
