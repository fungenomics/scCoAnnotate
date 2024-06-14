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
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import sys
# import pickle
import os
import random

# Set seed
random.seed(123456) 

#--------------- Parameters -------------------
ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
out_path = os.path.dirname(str(sys.argv[3]))
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
    
if 'batch' not in labels.columns:
  ## batch_key could not be None
  labels['batch'] = 'reference'
else:
  ## Just in case I add the prefix to be aware that is from the reference
  labels['batch'] = 'ref_' + labels['batch']
  print('@ Running with batches')


adata = ad.AnnData(X = ref,
                   obs = dict(obs_names=ref.index.astype(str),
                              labels=labels.label,
                              batch = labels.batch),
                   var = dict(var_names=ref.columns.astype(str))
                   )

## This step is necessary according to the authors
adata = remove_sparsity(adata)

batch_key = 'batch'
cell_type_key = 'labels'

#------------- Train scANVI -------------
# The data should be raw for scANVI and scVI
sca.models.SCVI.setup_anndata(adata,
                              batch_key=batch_key,
                              labels_key=cell_type_key)

vae = sca.models.SCVI(
    adata,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
#This step is slow, try to paralelize
# The 'train_dataloader' does not have many workers which may be a bottleneck. 
# Consider increasing the value of the `num_workers` argument` to `num_workers=9` in the `DataLoader` to improve performance.
vae.train()

scANVI_model = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
scANVI_model.train(max_epochs=20)

# Save the model to disk
print('@ SAVE MODEL')
scANVI_model.save(out_path,
                  overwrite=True)
print('@ DONE')
