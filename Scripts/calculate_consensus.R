
library(tidyverse)
library(glue)
args = commandArgs(trailingOnly = TRUE)
pred_path = args[1]
summary_path = args[2]
tools = strsplit(args[3], split = ' ')[[1]]
consensus_tools = strsplit(args[4], split = ' ')[[1]]
consensus_type = args[5]
ref_lab = args[6] 
print(args[7])
min_agree = as.numeric(strsplit(args[7], split = ' ')[[1]])
ontology_columns = args[8]
ontology_path = args[9]
f1_path = args[10]
CAWPE_type = strsplit(args[11], split = ' ')[[1]]
alpha = as.numeric(strsplit(args[12], split = ' ')[[1]])


print(tools)
print(consensus_tools)
print(min_agree)
print(ontology_path)
print(ontology_columns)
if(consensus_tools[1] == 'all'){
  consensus_tools = tools
}


print('CONSENSUS TOOLS: ')
print(consensus_tools)

print('CONSENSUS TYPE: ')
print(consensus_type)

## This function takes the data.frame ontology and a vector (pred)
## and convert the original label (from) to its ontology (to). 
## The output is the vector
apply_ontology <- function(df_ontology,
                           pred,
                           from = "labels",
                           to){
  ont <- setNames(nm = as.character(df_ontology[,from,drop=T]),
                  object = as.character(df_ontology[,to,drop=T])
  )
  match_ont <- as.character(ont[pred])
  pred[!is.na(match_ont)] <- match_ont[!is.na(match_ont)]
  return(pred)
}

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

files = list.files(pred_path, pattern = 'pred.csv', recursive = T, full.names = T)

l = list()
for(f in files){
   l[[basename(dirname(f))]] = data.table::fread(f) 
}

consensus = l %>% reduce(left_join, by = "cell") %>% rename('cellname' = 'cell') %>% select(all_of(c('cellname', tools))) %>% as.data.frame
rm(l)

if(consensus_type == 'majority' & all(min_agree != 0)){
  ref_labels = data.table::fread(ref_lab, header = T, fill=TRUE) %>% as.data.frame()
  if((ontology_columns != 'base') & (ontology_path != '')){
      ref_ontology = data.table::fread(ontology_path) %>% as.data.frame()
      ###Ontology step
      consensus[,tools] <- lapply(consensus[,tools,drop=F],
                          FUN = function(x){
                            apply_ontology(df_ontology = ref_ontology,
                                           pred = x,
                                           from = 'labels',
                                           to = ontology_columns)
                          }
      )
      print(consensus)
  }
  tmp = harmonize_unresolved(consensus, ref_labels)
  tmp = tmp %>% select(all_of(consensus_tools))
  for(mn_ag in min_agree){
    consensus[,paste0("Consensus_",as.character(mn_ag))] <- apply(tmp, 1, getmode, min = mn_ag)
  }
  rm(tmp)
}else if(consensus_type == "CAWPE" & all(alpha != 0)){
  ref_labels <- data.table::fread(glue('{dirname(ref_lab)}/labels_base.csv'), header = T, fill=TRUE) %>% as.data.frame()
  # get the paths to the probability files 
  prob_files = list.files(pred_path, 
                          pattern = paste0(paste0(consensus_tools, '_pred_score.csv'), collapse = '|'), 
                          recursive = T, 
                          full.names = T)
  
  for(CW_tp in CAWPE_type){
    if(CW_tp == 'CAWPE_T'){
      cols = c('tool')
    }else if(CW_tp == 'CAWPE_CT'){
      cols = c('tool', 'class')
    }
    # read F1 table, filter rows with consensus_tools and calculate mean F1
    F1 = read.csv(f1_path) %>%
      filter(tool %in% consensus_tools) %>%
      group_by(across(all_of(cols))) %>%
      summarize(F1 = mean(F1)) %>%
      ungroup()
    for(aph in alpha){
     # read probability matrixes
     data = lapply(prob_files, function(f){data.table::fread(f, check.names = F) %>% rename(cellname = V1) %>% mutate(tool = basename(dirname(f)))})
     
     # merge prob matrices and join with F1 scores
     data = bind_rows(data) %>% 
       select(cellname, all_of(ref_labels$label), tool) %>%
       pivot_longer(!c('cellname', 'tool'), 
                    names_to = 'class', 
                    values_to = 'prob') %>%
       left_join(F1, by = cols) %>%
       na.omit()
     
     data$CAWPE = apply(data, 1, CAWPE, alpha = aph)   
     # add up the CAWPE scores for each class for each cell 
     if(ontology_columns == 'base'){
      data = data %>%
         group_by(cellname, class) %>%
         summarise(CAWPE = sum(CAWPE)) %>%
         ungroup() %>% 
         group_by(cellname) %>% 
        dplyr::slice(which.max(CAWPE)) %>%
         rename(Consensus = class)  %>% 
         as.data.frame
      rownames(data) <- data$cellname
      # get final prediction by max CAWPE value 
      consensus[,paste0("CAWPE_",CW_tp,"_",aph)] = data[consensus$cellname,"CAWPE",drop=T]
      consensus[,paste0("Consensus_",CW_tp,"_",aph)] = data[consensus$cellname,"Consensus",drop=T]
     } else{
       ref_ontology = data.table::fread(ontology_path) %>% as.data.frame()
       print('@Here')
       ###Ontology step
       data = data %>%
         group_by(cellname, class) %>%
         summarise(CAWPE = sum(CAWPE)) %>%
         ungroup() %>% 
         as.data.frame
       data$ontology <- apply_ontology(df_ontology = ref_ontology,
                                       pred = data$class,
                                       from = 'labels',
                                       to = ontology_columns)
        data = data %>% 
         group_by(cellname, ontology) %>% 
         summarise(CAWPE = mean(CAWPE)) %>% 
         group_by(cellname) %>% 
         dplyr::slice(which.max(CAWPE)) %>%
         rename(Consensus = ontology)  %>% 
         as.data.frame
        rownames(data) <- data$cellname
        # get final prediction by max CAWPE value 
        consensus[,paste0("CAWPE_",CW_tp,"_",aph)] = data[consensus$cellname,"CAWPE",drop=T]
        consensus[,paste0("Consensus_",CW_tp,"_",aph)] = data[consensus$cellname,"Consensus",drop=T]
        consensus[,tools] <- lapply(consensus[,tools,drop=F],
                                    FUN = function(x){
                                      apply_ontology(df_ontology = ref_ontology,
                                                     pred = x,
                                                     from = 'labels',
                                                     to = ontology_columns)
                                    }
                                    )
     }
    }
  }
  
}

data.table::fwrite(consensus, summary_path, sep = '\t')

