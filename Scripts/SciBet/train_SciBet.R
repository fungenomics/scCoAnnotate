# load libraries and arguments 
library(data.table)
library(scibet)
library(WGCNA)
library(tidyverse)
library(Seurat)
library(glue)
set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
ref_path = args[1]
lab_path = args[2]
model_path = args[3]
threads = as.numeric(args[4])
out_path = dirname(model_path)
#--------------- Data -------------------

# read reference matrix and transpose 
message('@ READ REF')
### The matrix of the references is transpose since it needs to be normalize 
### with Seurat that expect a genes x cell.
ref <- data.table::fread(ref_path,
                         data.table=F,
                         header=T,
                         nThread=threads) %>%
  column_to_rownames('V1') %>% 
  WGCNA::transposeBigData()

message('@ DONE')

# read reference labels
labels <- data.table::fread(lab_path,
                            data.table=F,
                            header=T,
                            nThread=threads) %>%
  column_to_rownames('V1') 

# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(colnames(ref)))

# throw error if order is not the same 
if(!order){
  stop("@ Order of cells in reference and labels do not match")
}

### The matrix here is transposed since SciBet expect a cell x gene matrix and
### converted to data.frame since labels need to be add.
ref <- Seurat::NormalizeData(ref) %>% as.data.frame() %>% WGCNA::transposeBigData()

ref$label <- labels$label

#------------- Train SciBet -------------

Scibet_model <- scibet::Learn(expr = ref)

# save trained model 
message('@ SAVE MODEL')

save(Scibet_model,
     file = model_path)

message('@ DONE')

#------------- Other outputs --------------

### Transforming the model into a matrix
Scibet_matrix <- scibet::ExportModel(Scibet_model) %>% as.data.frame() %>% tibble::rownames_to_column(" ")

message('@ SAVE MODEL MATRIX')
data.table::fwrite(Scibet_matrix,
                   file = glue('{out_path}/model_matrix.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads
)
message('@ DONE')

### Taking the selected genes used in the model
selected_genes <- Scibet_matrix[,1,drop=T]
### Plotting in a marker heatmap
message('@ PLOTTING MARKER HEATMAP MATRIX')
pdf(glue('{out_path}/selected_markers_per_label.pdf',
         width = 12,
         height = 12)
)
marker_plot <- scibet::Marker_heatmap(expr = ref,
                                      gene = selected_genes)
plot(marker_plot)
dev.off()
message('@ DONE')

df_markers  <- marker_plot %>% .$data %>% filter(group == cell_type) %>% select(c('gene','cell_type','zscore'))
# write df_markers 
message('@ SAVE MARKERS')
data.table::fwrite(df_markers,
                   file = glue('{out_path}/selected_markers_per_label.csv'),
                   row.names = F,
                   col.names = T,
                   sep = ",",
                   nThread = threads
)
message('@ DONE')
#----------------------------------------