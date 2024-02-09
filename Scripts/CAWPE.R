
library(tidyverse)

#The arguments are provided in the appropriate snakefile
args = commandArgs(trailingOnly = TRUE)
pred_path = args[1]
summary_path = args[2]
tools = strsplit(args[3], split = ' ')[[1]]
consensus_tools = strsplit(args[4], split = ' ')[[1]]
consensus_type = args[5]
ref_lab = args[6]
f1_path = args[7]
alpha = as.numeric(args[8])

if(consensus_tools[1] == 'all'){
  consensus_tools = tools
}

if(consensus_type == 'CAWPE_T'){
  cols = c('tool')
}else if(consensus_type == 'CAWPE_CT'){
  cols = c('tool', 'class')
}

CAWPE = function(x, alpha = 4){
  (as.numeric(x['F1'])^alpha)*as.numeric(x['prob'])
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

print(prob_files)

# read probability matrixes
data = lapply(prob_files, function(f){data.table::fread(f, check.names = F) %>% rename(cell = V1) %>% mutate(tool = basename(dirname(f)))})

# read labels in matrix (used to filter the probbability matrix which sometimes contains additional information)
ref_labels = unique((data.table::fread(ref_lab, header = T, fill=TRUE))$label)

x = lapply(data, function(f){(f$cell)})

# 
data = bind_rows(data) %>% 
  select(cell, all_of(ref_labels), tool) %>%
  pivot_longer(!c('cell', 'tool'), 
               names_to = 'class', 
               values_to = 'prob') %>%
  left_join(F1, by = cols) %>%
  na.omit()

# calculate CAWPE score for each combination of cell, tool and class 
data$CAWPE = apply(data, 1, CAWPE, alpha = 4)

# get best prediction from each tool 
top_pred = data %>%
  group_by(cell, tool) %>%
  slice(which.max(prob)) %>%
  pivot_wider(id_cols = !c(prob, F1, CAWPE),
              names_from = tool, 
              values_from = class)

# add up the CAWPE scores for each class for each cell 
pred = data %>%
  group_by(cell, class) %>%
  summarise(CAWPE = sum(CAWPE)) %>%
  ungroup()

# get final prediction by max CAWPE value 
final_pred = pred %>% 
  group_by(cell) %>% 
  slice(which.max(CAWPE)) %>%
  rename(Consensus = class,
         cellname = cell) 

# save consensus 
data.table::fwrite(final_pred, summary_path, sep = '\t')

nrow(data)
