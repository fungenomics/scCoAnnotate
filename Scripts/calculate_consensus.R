
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
pred_path = args[1]
summary_path = args[2]
tools = strsplit(args[3], split = ' ')[[1]]
consensus_tools = strsplit(args[4], split = ' ')[[1]]
consensus_type = args[5]
ref_lab = args[6] 
f1_path = args[7]
alpha = as.numeric(args[8])

print(tools)
print(consensus_tools)

if(consensus_tools[1] == 'all'){
  consensus_tools = tools
}

if(consensus_type == 'majority'){
  min_agree = 1 
}else{
  min_agree = as.numeric(consensus_type) 
}

print('CONSENSUS TOOLS: ')
print(consensus_tools)

print('CONSENSUS TYPE: ')
print(consensus_type)

harmonize_unresolved = function(pred, ref_labels){
  pred %>%
  column_to_rownames('cellname') %>%
  mutate(across(where(is.character), ~ifelse(. %in% c(ref_labels$label), ., 'Unresolved'))) %>%
  rownames_to_column('cellname') %>%
  return()
}

getmode = function(v, min = 1){
  uniqv = unique(v)
  matches = tabulate(match(v, uniqv))
  max_match = max(matches)
  
  ties = ifelse(length(which(matches == max_match)) > 1, T, F)
  
  if (max_match < min) {
    return("No Consensus")
  } else if (ties) {
    return("No Consensus") 
  } else {
    return(uniqv[which.max(matches)])
  }
}

CAWPE = function(x, alpha = 4){
  (as.numeric(x['F1'])^alpha)*as.numeric(x['prob'])
}

ref_labels = data.table::fread(ref_lab, header = T, fill=TRUE) %>% column_to_rownames('V1')

files = list.files(pred_path, pattern = 'pred.csv', recursive = T, full.names = T)

l = list()
for(f in files){
   l[[basename(dirname(f))]] = data.table::fread(f) 
}

consensus = l %>% reduce(left_join, by = "cell") %>% rename('cellname' = 'cell') %>% select(all_of(c('cellname', tools)))
rm(l)

if(consensus_type == 'majority'){

  tmp = harmonize_unresolved(consensus, ref_labels)
  tmp = tmp %>% select(all_of(consensus_tools))
  consensus$Consensus = apply(tmp, 1, getmode, min = min_agree)
  rm(tmp)

}else if(consensus_type %in% c('CAWPE_T', 'CAWPE_CT')){
  
  if(consensus_type == 'CAWPE_T'){
     cols = c('tool')
  }else if(consensus_type == 'CAWPE_CT'){
    cols = c('tool', 'class')
  }
  
  # read F1 table, filter rows with consensus_tools and calculate mean F1
  F1 = read.csv(f1_path) %>%
       filter(tool %in% consensus_tools) %>%
       group_by(across(all_of(cols))) %>%
       summarize(F1 = mean(F1)) %>%
       ungroup()

  # get the paths to the probability files 
  prob_files = list.files(pred_path, 
                        pattern = paste0(paste0(consensus_tools, '_pred_score.csv'), collapse = '|'), 
                        recursive = T, 
                        full.names = T)

  # read probability matrixes
  data = lapply(prob_files, function(f){data.table::fread(f, check.names = F) %>% rename(cellname = V1) %>% mutate(tool = basename(dirname(f)))})

  # read labels in matrix (used to filter the probbability matrix which sometimes contains additional information)
  ref_labels = unique((data.table::fread(ref_lab, header = T, fill=TRUE))$label)

  # merge prob matrices and join with F1 scores
  data = bind_rows(data) %>% 
         select(cellname, all_of(ref_labels), tool) %>%
         pivot_longer(!c('cellname', 'tool'), 
                      names_to = 'class', 
                      values_to = 'prob') %>%
         left_join(F1, by = cols) %>%
         na.omit()

  data$CAWPE = apply(data, 1, CAWPE, alpha = 4)

  # add up the CAWPE scores for each class for each cell 
  data = data %>%
         group_by(cellname, class) %>%
         summarise(CAWPE = sum(CAWPE)) %>%
         ungroup() %>% 
         group_by(cellname) %>% 
         slice(which.max(CAWPE)) %>%
         rename(Consensus = class) 

   # get final prediction by max CAWPE value 
  consensus = consensus %>% 
              left_join(data, by = 'cellname') 

}

data.table::fwrite(consensus, summary_path, sep = '\t')

