# load libraries and arguments 

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
model_path = args[3]
threads = as.numeric(args[4])

# path for other outputs (depends on tools)
out_path = dirname(model_path)

#--------------- Data -------------------

# read reference matrix 
message('@ READ REF')
ref = data.table::fread(ref_path, nThread=threads, header=T, data.table=F) %>%
      column_to_rownames('V1')
message('@ DONE')

# read reference labels
labels = data.table::fread(lab_path, header=T, data.table=F) %>%
         column_to_rownames('V1')

# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(rownames(ref)))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# CODE FOR ANY DATA MANIPULATION NEEDED BEFORE TRAINING HERE 


#------------- Train <tool_name> -------------

# CODE FOR TRAINING HERE  

# save trained model 
message('@ SAVE MODEL')
save(<tool_model>, file = model_path)
message('@ DONE')

#---------------------------------------------
