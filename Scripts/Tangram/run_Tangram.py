#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 17:11:20 2023

@author: tomas.vega
"""
## Python 3.11.2
#--------------- Libraries -------------------
import tangram as tg
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import sys
import random
### Set seed
random.seed(123456) 

#--------------- Parameters -------------------
ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
sp_dataset = str(sys.argv[4])
mode = str(sys.argv[5])
#if str(sys.argv[5]) != 'None':
#  cluster_label = str(sys.argv[5])
#else:
#  cluster_label = None

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

ad_sc = ad.AnnData(X = ref,
                   obs = dict(obs_names=ref.index.astype(str),
                              label = labels['label']),
                   var = dict(var_names=ref.columns.astype(str))
                   )

## Now I normalize the matrix with scanpy:
#Normalize each cell by total counts over all genes,
#so that every cell has the same total count after normalization.
#If choosing `target_sum=1e6`, this is CPM normalization
#1e4 similar as Seurat
sc.pp.normalize_total(ad_sc, target_sum=1e4)
#Logarithmize the data:
sc.pp.log1p(ad_sc)

# read spatial query
ad_sp = sc.read_h5ad(sp_dataset)

### Use only the same genes in the spatial and the query
tg.pp_adatas(ad_sc, ad_sp, genes=None)

if mode == 'clusters':
    ad_map = tg.map_cells_to_space(
                   ad_sc,
                   ad_sp,
                   mode=mode,
                   cluster_label='label')
else:
    ad_map = tg.map_cells_to_space(
                   ad_sc,
                   ad_sp,
                   mode=mode)


tg.project_cell_annotations(ad_map, ad_sp, annotation='label')
df = ad_sp.obsm['tangram_ct_pred'].copy()

# Function to get the column name and maximum value for each row
def get_max_column_and_value(row):
    pred_label = row.idxmax()
    proba_label = row.max()
    return pred_label, proba_label

# Create a new column 'max_column' with the column name containing the maximum value for each row
df['pred_label'], df['score_label'] = zip(*df.apply(get_max_column_and_value, axis=1))
# Create a new column 'unknown_max_column' to store 'max_column' as 'unknown' if 'max_value' is lower than the threshold
#df['pred_label_reject'] = df.apply(lambda row: 'Unknown' if row['proba_label'] < threshold else row['pred_label'], axis=1)

print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'spot': df.index, 'Tangram': df.pred_label})
pred_df.to_csv(out_path, index = False)
print('@ DONE')

#------------- Other outputs --------------
### Save transfered score matrix
print('@ WRITTING TRANSFERED SCORE MATRIX ')
filename = out_other_path + '/label_score_matrix.csv'
df.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')

#### Save map object that contains the cell x spot prob
print('@ WRITTING OUTPUT OBJECT ')
filename = out_other_path + '/Tangram_mapped_object.h5ad'
ad_map.write_h5ad(filename= filename,
                  compression='gzip')
print('@ DONE ')

#### Save cell x spot prob matrix
print('@ WRITTING PROB MATRIX ')
filename = out_other_path + '/prob_matrix.h5ad'
prob_matrix = pd.DataFrame(ad_map.X,index = ad_map.obs_names,columns = ad_map.var_names)
prob_matrix.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')
