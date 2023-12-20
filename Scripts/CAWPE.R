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

if(consensus_tools[1] == 'all'){
  consensus_tools = tools
}

CAWPE = function(v, cva, ref_labels, consensus_tools, column_names, alpha = 4){
  #Slash is added to create exact match for labels
  tmp_labels = paste0(ref_labels,"/")
  
  CAWPE_values = list()
  
  #for each label 
  for (lab in tmp_labels){
    
    # get the F1 scores for each tool for one lable 
    cva_one_label = cva %>% filter(class == lab)
    
    # get the column names containing the label (tool.label format)?? 
    # can be done with a paste instead maybe 
    # one_label_prob_columns = paste0(consensus_tools, '.', lab)
    one_label_prob_columns = grep(lab, column_names)

    # filtered for one label (named named vector)
    v_one_label = as.list(v)[one_label_prob_columns]
    
    CAWPE_values[lab] = 0
    
    for(c_tool in consensus_tools){
      #There is now only one column for the label-tool pair
      label_tool_prob = as.double(v_one_label[grep(c_tool, names(v_one_label))])
      
      # check if there is a value 
      if(length(label_tool_prob) == 0){label_tool_prob = 0}
      
      # calculate value per tool and lab
      tool_contribution = (cva_one_label %>% filter(tool == c_tool))[,"F1"]^alpha*label_tool_prob
      
      # add to sum 
      CAWPE_values[lab] = CAWPE_values[lab] + tool_contribution
    }
    
  }
  #Ref_labels iteration determines the order anyhow, this extracts the label with no slash
  return(ref_labels[which.max(unlist(CAWPE_values))])
  
}


F1_file = read.csv(f1_path)

#Process the F1 data
#Add slash to distinguish labels
F1_file$class = paste0(F1_file$class,"/")

#Keep only rows with active consensus_tools
F1_filtered = F1_file %>%
  filter(tool %in% consensus_tools)

#Calculate the mean F1 score for each class and tool (this is the cross-validation accuracy estimate)
cva = F1_filtered %>%
  group_by(tool, class) %>%
  summarize(F1 = mean(F1))


#The filtering assumes that tools have a parent directory with the consensus tool name

unfiltered_pred_files =list.files(pred_path, pattern = 'pred.csv', recursive = T, full.names = T)
pred_files = grep(paste(paste0("/",consensus_tools,"/"),collapse="|"), unfiltered_pred_files, value=TRUE)
rm(unfiltered_pred_files)

unfiltered_prob_files =list.files(pred_path, pattern = 'prob_matrix.csv', recursive = T, full.names = T)
prob_files = grep(paste(paste0("/",consensus_tools,"/"),collapse="|"), unfiltered_prob_files, value=TRUE)
rm(unfiltered_prob_files)

#If we do not use unique, we get all the training data redundancies
ref_labels = unique((data.table::fread(ref_lab, header = T))$label)

l = list()
l2 = list()
for(f in prob_files){
  l[[basename(dirname(f))]] = data.table::fread(f) 
  }
for(f in pred_files){
  l2[[basename(dirname(f))]] = data.table::fread(f) 
}

#cbind automatically provides the tool name to the table columns! e.g. scPred.ASEP
consensus_calc_table = do.call(cbind,l)

#Add slash to distinguish between labels that are similar
colnames(consensus_calc_table) = paste0(colnames(consensus_calc_table),"/")
column_names = colnames(consensus_calc_table)
consensus = l2 %>% reduce(left_join, by = "cell") %>% rename('cellname' = 'cell')
rm(l)
rm(l2)

#Append the consensus prediction to the table of all best predictions
consensus$Consensus = apply(consensus_calc_table, 1, CAWPE, cva, ref_labels, consensus_tools, column_names)
data.table::fwrite(consensus, summary_path, sep = '\t')
