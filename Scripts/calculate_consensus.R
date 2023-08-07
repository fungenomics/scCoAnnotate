
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
pred_path = args[1]
summary_path = args[2]
tools = strsplit(args[3], split = ' ')[[1]]
consensus_tools = strsplit(args[4], split = ' ')[[1]]
ref_lab = args[5] 

print(tools)
print(consensus_tools)

if(consensus_tools[1] == 'all'){
  consensus_tools = tools
}

print(consensus_tools)

harmonize_unsure = function(pred, ref_labels){
  pred %>%
  column_to_rownames('cellname') %>%
  mutate(across(where(is.character), ~ifelse(. %in% c(ref_labels$label), ., 'Unsure'))) %>%
  rownames_to_column('cellname') %>%
  return()
}

getmode = function(v) {
  uniqv = unique(v)
  matches = tabulate(match(v, uniqv))
  max_match = max(matches)
  
  ties = ifelse(length(which(matches == max_match)) > 1, T, F)
  
  if (max_match == 1) {
    return("No Consensus")
  } else if (ties) {
    return("No Consensus") 
  } else {
    return(uniqv[which.max(matches)])
  }
}

ref_labels = data.table::fread(ref_lab, header = T) %>% column_to_rownames('V1')

files =list.files(pred_path, pattern = 'pred.csv', recursive = T, full.names = T)

l = list()
for(f in files){
   l[[basename(dirname(f))]] = data.table::fread(f) 
}

consensus = l %>% reduce(left_join, by = "cell") %>% rename('cellname' = 'cell')
rm(l)

tmp = harmonize_unsure(consensus, ref_labels)
tmp = tmp %>% select(all_of(consensus_tools))
consensus$Consensus = apply(tmp, 1, getmode)
rm(tmp)

data.table::fwrite(consensus, summary_path, sep = '\t')


