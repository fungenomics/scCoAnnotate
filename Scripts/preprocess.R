
library(tidyverse)

initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]) 

source(paste0(dirname(script.name), "/Functions/functions.R"))

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
query_paths = strsplit(args[2], split = ' ')[[1]]
out = args[3]
convert_genes = as.logical(args[4])
lab_path = args[5]
reference_name = args[6]
query_names = strsplit(args[7], split = ' ')[[1]]

min_cells = as.numeric(args[8])
if(is.na(min_cells)){
  stop("The minimum number of cells specified is not a numeric value")
}

downsample_value = as.numeric(args[9]) 
if(is.na(downsample_value)){
  stop("The downsample value specified is not a numeric value")
}

downsample_per_class = as.logical(args[10])
if(is.na(downsample_per_class)){
  stop("The downsample stratified specified is not a logical value")
}

names(query_paths) = query_names

batch_path = args[11]
if(batch_path == 'None'){
  batch_path = NULL
}

print(batch_path)
# ----- PREPROCESS REFERENCE ----------------------
tmp <- get_data_reference(ref_path = ref_path,
                          lab_path = lab_path,
                          batch_path = batch_path)
data <- list()
data[['ref']] <- tmp$exp
lab           <- tmp$lab
rm(tmp)

# downsample 
if(downsample_value != 0){
  lab = downsample(lab, downsample_per_class, downsample_value)
}

# remove small clusters 
if(min_cells > 0){
  lab = remove_small_clusters(lab, min_cells)
}
# filter reference for donwsampled cells 
data[['ref']] = data[['ref']][rownames(lab),]

# save downsampled lables 
save.df <- data.frame(cells= rownames(lab), 
                      lab)

colnames(save.df)[1] <- ""

data.table::fwrite(save.df,
                   file = paste0(out, '/model/', reference_name, '/downsampled_labels.csv'),
                   col.names = T,
                   row.names=F,
                   sep = ",")
rm(save.df)

# if specified by user, convert reference gene names from mouse to human
if(convert_genes){
  
  message('@ CONVERTING GENE NAMES')

  # include functions and libraries for conversion
  library(Orthology.eg.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(WGCNA)

  # convert
  hg = mapfun(colnames(data[['ref']])) %>% dplyr::select(Mouse_symbol, Human_symbol) 
  
  # output list of mouse genes that were not converted
  not_converted = hg %>% filter(is.na(Human_symbol)) %>% .$Mouse_symbol
  
  data.table::fwrite(as.list(not_converted), file = paste0(out, '/model/', reference_name, '/genes_not_converted.csv'), sep = ',')

  # throw error if more than threshold % genes not converted
  threshold = 0.5
  if(length(not_converted) > threshold*length(colnames(data[['ref']]))){
    stop(paste0("@ More than ",threshold*100,"% of mouse genes in reference could not be converted to human"))
  }

  # modify reference matrix to contain converted genes
  data[['ref']] = data[['ref']] %>%
    transposeBigData() %>%
    rownames_to_column('Mouse_symbol') %>%
    inner_join(hg %>% filter(!is.na(Human_symbol)), 
               by = 'Mouse_symbol') %>%
    dplyr::select(-Mouse_symbol) %>%
    column_to_rownames('Human_symbol') %>%
    transposeBigData() 

}

# ----- QUERY --------------------------------

# read query 
for(i in 1:length(query_paths)){
  
  print(query_paths[i])
  tmp = get_data_query(query_path = query_paths[i])
  query = names(query_paths)[i]
  
  print(query)
  data[[query]] = tmp
}

# ----- GENE INTERSECT ------------------------

# get genes for each data frame (colnames)
genes = lapply(data, function(x){(colnames(x))})

# reduce set of genes to the intersect 
common_genes = Reduce(intersect,genes)
print(paste0('@Found ', length(common_genes), ' in common'))

# throw error if number of common genes below % threshold of genes in any of provided datasets (ref or query) 
threshold = 0.25
frac = lapply(genes, function(x){length(common_genes)/length(x)})

if(any(frac < threshold)){
  names(frac) = names(data)
  print(frac)
  stop(paste0("@ In at least one provided dataset (ref or query), less than ",threshold*100,"% of genes appear in common gene set. See above for the fraction of genes from each dataset appearing in common gene set (note: samples with few genes will have higher fractions)"))
}

# save common genes 
data.table::fwrite(data.frame('common_genes' = common_genes), file = paste0(out, '/model/', reference_name, '/common_genes.csv'))

# filter each data set for common genes
data = lapply(data, function(x){x[,common_genes]})

#----- SAVE DATA ----------------------------------------

# save reference 
tmp = data[['ref']] %>% rownames_to_column()
colnames(tmp)[1] = " "

data.table::fwrite(tmp, file = paste0(out, '/model/', reference_name, '/expression.csv'), sep = ',')

# save query 
query_names = names(data)[!names(data) == 'ref']
for(q in query_names){
  print(q)
  
  tmp = data[[q]] %>% rownames_to_column()
  colnames(tmp)[1] = " "
  
  data.table::fwrite(tmp, file = paste0(out, '/', q, '/', reference_name, '/expression.csv'), sep = ',')
}

#---------------------------------------------------------

