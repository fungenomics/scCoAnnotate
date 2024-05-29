#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28

@author: tomas.vega
"""
## Python 3.11.2
#--------------- Libraries -------------------
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
from scarches.models.scpoli import scPoli
import sys
import pickle
import os
import random

# Set seed
random.seed(123456) 

#--------------- Parameters -------------------
ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
out_path = str(sys.argv[3])

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
    
adata = ad.AnnData(X = ref,
                   obs = dict(obs_names=ref.index.astype(str),
                              labels=labels.label),
                   var = dict(var_names=ref.columns.astype(str))
                   )

# Now I normalize the matrix with scanpy:
# Normalize each cell by total counts over all genes,
# so that every cell has the same total count after normalization.
# If choosing `target_sum=1e6`, this is CPM normalization
# 1e4 similar as Seurat
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data:
sc.pp.log1p(adata)

## condition_key could not be None
adata.obs['condition'] = 'reference'

condition_key = 'condition'
cell_type_key = 'labels'
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
#------------- Train scPoli -------------
scpoli_model = scPoli(
    adata=adata,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    embedding_dims=5,
    recon_loss='nb',
)

scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)
# Save the model to disk
print('@ SAVE MODEL')
pickle.dump(scpoli_model, open(out_path, 'wb'))
print('@ DONE')
