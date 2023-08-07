
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
query_paths = args[2]
query_paths = strsplit(query_paths, split = ' ')[[1]]
out_path = args[3]

paths = c(ref_path, query_paths)
print(paths)

l = list()
for(i in 1:length(paths)){
  tmp = data.table::fread(paths[i], header = T) %>% column_to_rownames('V1')
  l[[i]] = tmp
}

# get genes for each data frame (colnames)
genes = lapply(l, function(x){(colnames(x))})

# reduce set of genes to the intersect 
genes = Reduce(intersect, genes)

# save common genes 
data.table::fwrite(data.frame('common_genes' = genes))

# filter expression data to the intersect 
for(i in 1:length(l)){
  tmp = l[[i]][,genes]
  
  # calculate fraction of genes kept 
  fr = (ncol(tmp)/ncol(l[[i]]))

  # throw error if fraction is less than 0.5  
  if(fr < 0.5){
    stop(paste(" ", 
    	"@ Number of genes after filtering is less than 50% of original data for sample:", 
    	paths[i], 
    	sep="\n"))
  }

  # save filtered data 
  if(i == 1){
     data.table::fwrite(tmp, )
  }else{
    data.table::fwrite(tmp)
  }
}
