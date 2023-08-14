
library(tidyverse)
library(WGCNA)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
query_paths = strsplit(args[2], split = ' ')[[1]]
out_path = args[3]
convert_genes = as.logical(args[4])

l = list()

# read reference 
l[['ref']] = data.table::fread(ref_path, header = T) %>% column_to_rownames('V1')

# read query 
for(p in query_paths){
  tmp = data.table::fread(p, header = T) %>% column_to_rownames('V1')
  query = basename(dirname(p))
  l[[query]] = tmp
}

# if specified by user, convert reference gene names from mouse to human
if(convert_genes){

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
  data.table::fwrite(as.list(not_converted), file = paste0(out_path, '/genes_not_converted.csv'), sep = ',')

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
