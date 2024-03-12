
library(tidyverse)
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
downsample_stratified = as.logical(args[10])
if(is.na(downsample_stratified)){
  stop("The downsample stratified specified is not a logical value")
}
downsample_stratified = if(downsample_stratified) "label" else NULL

names(query_paths) = query_names

l = list()

# read reference 
l[['ref']] = data.table::fread(ref_path, header = T) %>% column_to_rownames('V1')
#read labels 
lab = data.table::fread(lab_path, header = T)
#If the value is > 1 it will use it as an absolute value
#If the value is a fraction it will downsample the fraction.
# If the group has lower number of cells than the downsample_value the total 
#amount of cells are retain for that cluster
if(downsample_value != 0){
  if(downsample_value >= 1){
      lab = lab %>%  group_by(across(all_of(downsample_stratified))) %>% 
        dplyr::slice_sample(n = downsample_value,replace = F)
  } else{
      lab = lab %>% group_by(across(all_of(downsample_stratified))) %>% 
        dplyr::slice_sample(prop = downsample_value,replace = F) 
  }
  
  l[['ref']] = l[['ref']][lab$V1,]
}
#do the convertion to rownames after since the column is needed to downsampling
lab = lab %>% column_to_rownames('V1')

if(min_cells > 0){
  rmv_labels = names(which(table(lab$label) < min_cells))
  lab = lab %>% filter((!label %in% rmv_labels))
  message(paste0(paste0(rmv_labels,collapse = '-'),' classes were remove because of lower number of cells (< ',as.character(min_cells),')'))
  #filtering the cells from the filtered classes
  l[['ref']] = l[['ref']][rownames(lab),]
}

save.df <- data.frame(cells= rownames(lab),
                      label= lab$label)
colnames(save.df)[1] <- ""
data.table::fwrite(save.df,
                   file = paste0(out, '/model/', reference_name, '/downsampled_labels.csv'),
                   col.names = T,
                   row.names=F,
                   sep = ",")
rm(save.df)

# read query 
for(i in 1:length(query_paths)){
  print(query_paths[i])
  tmp = data.table::fread(query_paths[i], header = T) %>% column_to_rownames('V1')
  query = names(query_paths)[i]
  print(query)
  l[[query]] = tmp
}

# if specified by user, convert reference gene names from mouse to human
if(convert_genes){
  
  message('@ CONVERTING GENE NAMES')

  # include functions and libraries for conversion
  library(Orthology.eg.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(WGCNA)

  mapfun = function(mousegenes){
    gns    = mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
    mapped = AnnotationDbi::select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
    naind  = is.na(mapped$Homo_sapiens)
    hsymb  = mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
    out    = data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
    out$Human_symbol[!naind] = hsymb
    return(out)
  }

  # convert
  hg = mapfun(colnames(l[['ref']])) %>% dplyr::select(Mouse_symbol, Human_symbol) 
  
  # output list of mouse genes that were not converted
  not_converted = hg %>% filter(is.na(Human_symbol)) %>% .$Mouse_symbol
  data.table::fwrite(as.list(not_converted), file = paste0(out, '/model/', reference_name, '/genes_not_converted.csv'), sep = ',')

  # throw error if more than threshold % genes not converted
  threshold = 0.5
  if(length(not_converted) > threshold*length(colnames(l[['ref']]))){
    stop(paste0("@ More than ",threshold*100,"% of mouse genes in reference could not be converted to human"))
  }

  # modify reference matrix to contain converted genes
  l[['ref']] = l[['ref']] %>%
    transposeBigData() %>%
    rownames_to_column('Mouse_symbol') %>%
    inner_join(hg %>% filter(!is.na(Human_symbol)), 
               by = 'Mouse_symbol') %>%
    dplyr::select(-Mouse_symbol) %>%
    column_to_rownames('Human_symbol') %>%
    transposeBigData() 

}

# get genes for each data frame (colnames)
genes = lapply(l, function(x){(colnames(x))})

# reduce set of genes to the intersect 
common_genes = Reduce(intersect,genes)
print(paste0('@Found ', length(common_genes), ' in common'))

# throw error if number of common genes below % threshold of genes in any of provided datasets (ref or query) 
threshold = 0.25
frac = lapply(genes, function(x){length(common_genes)/length(x)})

if(any(frac < threshold)){
  names(frac) = names(l)
  print(frac)
  stop(paste0("@ In at least one provided dataset (ref or query), less than ",threshold*100,"% of genes appear in common gene set. See above for the fraction of genes from each dataset appearing in common gene set (note: samples with few genes will have higher fractions)"))
}

# save common genes 
data.table::fwrite(data.frame('common_genes' = common_genes), file = paste0(out, '/model/', reference_name, '/common_genes.csv'))

# filter each data set for common genes
l = lapply(l, function(x){x[,common_genes]})

# save reference 
tmp = l[['ref']] %>% rownames_to_column()
colnames(tmp)[1] = " "
data.table::fwrite(tmp, file = paste0(out, '/model/', reference_name, '/expression.csv'), sep = ',')

# save query 
query_names = names(l)[!names(l) == 'ref']
for(q in query_names){
  print(q)
  tmp = l[[q]] %>% rownames_to_column()
  colnames(tmp)[1] = " "
  data.table::fwrite(tmp, file = paste0(out, '/', q, '/', reference_name, '/expression.csv'), sep = ',')
}

# save unique labels (for downstream report color pal)
# lab = data.table::fread(lab_path, header = T) %>% column_to_rownames('V1')
lab = data.frame(label = unique(lab$label))
data.table::fwrite(lab, file = paste0(out, '/model/', reference_name, '/labels.csv'), sep = ',')

