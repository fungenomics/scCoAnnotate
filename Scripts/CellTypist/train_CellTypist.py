## Python 3.11.2
#--------------- Libraries -------------------
import numpy as np
import pandas as pd
import celltypist
import anndata as ad
import sys
import scanpy as sc
import pickle
import os
import random

### Set seed
random.seed(123456) 

#--------------- Parameters -------------------

ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threads = int(sys.argv[4])
feature_selection = bool(sys.argv[5])

#--------------- Data -------------------------

# read the data
ref = pd.read_csv(ref_path,
                  index_col=0,
                  sep=',',
                  engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

labels = pd.read_csv(lab_path,
                     index_col = 0,
                     sep=',')

# check if cell names are in the same order in labels and ref
order = all(labels.index == ref.index)

# throw error if order is not the same 
if not order:
  sys.exit("@ Order of cells in reference and labels do not match")

# create AnnData object 
adata = ad.AnnData(X = ref,
                  obs = dict(obs_names=ref.index.astype(str),
                             label = labels['label']),
                  var = dict(var_names=ref.columns.astype(str))
                  )

# normalize the matrix with scanpy:
# normalize each cell by total counts over all genes,
# so that every cell has the same total count after normalization.
# If choosing `target_sum=1e6`, this is CPM normalization
# 1e4 similar as Seurat
sc.pp.normalize_total(adata, target_sum=1e4)

#Logarithmize the data:
sc.pp.log1p(adata)

#------------- Train CellTypist -------------

model = celltypist.train(adata,
                         labels = "label",
                         transpose_input = False,
                         check_expression = True,
                         feature_selection = feature_selection,
                         n_jobs = threads)


#Save the model to disk
print('@ SAVE MODEL')                             
model.write(out_path)
print('@ DONE')

#------------- Other outputs --------------

# We can extract the top markers, I get the top 10 for each cell-type applying the
# function extract_top_markers
dataframes = []
for cell_type in model.cell_types:
    top_markers = model.extract_top_markers(cell_type, 10)
    
    # Create a DataFrame for the current cell type's top markers
    df = pd.DataFrame(top_markers, columns=['Marker'])
    df['CellType'] = cell_type  # Add a column to store the cell type
    
    # Append the DataFrame to the list
    dataframes.append(df)

# Concatenate all DataFrames into a single DataFrame
markers_df = pd.concat(dataframes, ignore_index=True)

filename = out_other_path + "/top10_model_markers_per_celltype.csv"
markers_df.to_csv(filename,
                  index=False) 
