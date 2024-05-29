#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28

@author: tomas.vega
"""
## Python 3.11.2
#--------------- Libraries -------------------
import pandas as pd
from scarches.models.scpoli import scPoli
import anndata as ad
import sys
import os
import scanpy as sc
import pickle
import random
import numpy as np
# Set seed
random.seed(123456) 

# Function to get the column name and maximum value for each row
def get_max_column_and_value(row):
    pred_label = row.idxmax()
    proba_label = row.max()
    return pred_label, proba_label
  
# def normalize(x):
#   """Compute softmax values for each sets of scores in x."""
#   return (x - np.min(x) / np.range(x))

def normalize(x):
    """Normalize each row to the range 0-1."""
    min_val = np.min(x, axis=1, keepdims=True)
    max_val = np.max(x, axis=1, keepdims=True)
    return (x - min_val) / (max_val - min_val)
  
#--------------- Parameters -------------------
sample_path = str(sys.argv[1])
model_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))

#--------------- Data -------------------

# read query matrix
print('@ READ QUERY')
query = pd.read_csv(sample_path,
                    index_col=0,
                    sep=',',
                    engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

print('@ DONE')

# load model 
print('@ LOAD MODEL')
scPoli_model = pickle.load(open(model_path, 'rb'))
print('@ DONE')

# Query preprocessing
query = ad.AnnData(X = query,
                   obs = dict(obs_names=query.index.astype(str)),
                   var = dict(var_names=query.columns.astype(str))
                   )

# Now I normalize the matrix with scanpy:
# Normalize each cell by total counts over all genes,
# so that every cell has the same total count after normalization.
# If choosing `target_sum=1e6`, this is CPM normalization
# 1e4 similar as Seurat
sc.pp.normalize_total(query, target_sum=1e4)

# Logarithmize the data:
sc.pp.log1p(query)

query.obs['condition'] = 'query'
query.obs['labels'] = 'Unknown'

#----------- Predict scPoli --------
## Need the intersection
scpoli_query = scPoli.load_query_data(
    adata=query,
    reference_model=scPoli_model,
    labeled_indices=[],
)

scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)

pred_proba = scpoli_query.classify(query,
                                   scale_uncertainties=False,
                                   log_distance=False)
                                   
print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'cell': query.obs_names,
                       'scPoli': pred_proba['labels']['preds']})
                       
pred_df.to_csv(out_path,
               index = False)
print('@ DONE')


#----------------------------------------
# Get the distances and reescale to an score that the higher the better
mtx = 1/(pred_proba['labels']['weighted_distances'] + 1)
mtx = normalize(mtx)
# mtx = pred_proba['labels']['weighted_distances']
df = pd.DataFrame(mtx,index = query.obs_names,columns = scPoli_model.model_cell_types)
df['pred_label'], df['score_label'] = zip(*df.apply(get_max_column_and_value, axis=1))
all(pred_proba['labels']['preds'] == df.pred_label)
# Save the score matrix
print('@ WRITTING PROB MATRIX ')
filename = out_other_path + '/scPoli_pred_score.csv'
df.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')
