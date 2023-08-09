
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
query_paths = strsplit(args[2], split = ' ')[[1]]
out_path = args[3]

l = list()

# read reference 
l[['ref']] = data.table::fread(ref_path, header = T) %>% column_to_rownames('V1')

# read query 
for(p in query_paths){
  tmp = data.table::fread(p, header = T) %>% column_to_rownames('V1')
  query = basename(dirname(p))
  l[[query]] = tmp
}

# get genes for each data frame (colnames)
genes = lapply(l, function(x){(colnames(x))})

# reduce set of genes to the intersect 
genes = Reduce(intersect, genes)

# save common genes 
data.table::fwrite(data.frame('common_genes' = genes), paste0(out_path, '/common_genes.csv'))

# filter each data set for common genes
l = lapply(l, function(x){x[,genes]})

# save reference 
tmp = l[['ref']] %>% rownames_to_column()
colnames(tmp)[1] = " "
data.table::fwrite(tmp, file = paste0(out_path, '/expression.csv'), sep = ',')

# save query 
query_names = names(l)[!names(l) == 'ref']
for(q in query_names){
  print(q)
  tmp = l[[q]] %>% rownames_to_column()
  colnames(tmp)[1] = " "
  data.table::fwrite(tmp, file = paste0(out_path, '/', q, '/expression.csv'), sep = ',')
}
