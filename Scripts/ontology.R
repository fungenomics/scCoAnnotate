
library(tidyverse)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)

out = args[1]
reference_name = args[2]
lab_path = args[3]
ontology_path = args[4]                                                                                                                            
ontology_columns = strsplit(args[5], split = ' ')[[1]]

print(out)
print(lab_path)
print(reference_name)
print(ontology_path)
print(ontology_columns)

#----- SAVE ONTOLOGY -----------------------------------

dir.create(paste0(out, '/model/', reference_name, '/ontology/'), recursive = T)

if(length(ontology_columns) == 1 & ontology_columns[1] == 'label'){
  
  lab = data.table::fread(lab_path, header = T)
  print(lab)

  ont = data.frame(label = unique(lab$label))
  print(ont)

  data.table::fwrite(ont,
                     file = paste0(out, '/model/', reference_name, '/ontology/ontology.csv'),
                     sep = ',')
}else{
  ont = data.table::fread(ontology_path) 
  print(ont)

  data.table::fwrite(ont,
                     file = paste0(out, '/model/', reference_name, '/ontology/ontology.csv'),
                     sep = ',')
}

#--------------------------------------------------------
