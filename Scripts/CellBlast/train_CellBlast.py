#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 2024

@author: tomas.vega
"""
## Python 3.11.2
#--------------- Libraries -------------------
import pandas as pd
import anndata as ad
import Cell_BLAST as cb
import sys
import os
import random

# Set seed
random.seed(123456) 

#--------------- Parameters -------------------
ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
out_path = os.path.dirname(str(sys.argv[3]))
threads = int(sys.argv[4])
n_models = int(sys.argv[5])
cb.config.N_JOBS = threads
cb.config.RANDOM_SEED = 123456

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

batch = None
if 'batch' not in labels.columns:
  ## batch_key could not be None
  labels['batch'] = 'reference'
else:
  ## Just in case I add the prefix to be aware that is from the reference
  labels['batch'] = 'ref_' + labels['batch']
  batch = 'batch'
  print('@ Running with batches')
  
adata = ad.AnnData(X = ref,
                   obs = dict(obs_names=ref.index.astype(str),
                              labels=labels.label,
                              batch = labels.batch),
                   var = dict(var_names=ref.columns.astype(str))
                   )


cb.data.find_variable_genes(adata,
                            grouping=batch)

#------------- Train CellBlast -------------
print('@ Training CellBlast')
models = []
for i in range(n_models):
  print(i)
  models.append(cb.directi.fit_DIRECTi(
                adata,
                genes=adata.var.query("variable_genes").index,
                latent_dim=10,
                cat_dim=20,
                random_seed=i
    ))

blast_model = cb.blast.BLAST(models, adata)
print('@ DONE')

# Save the model to disk
print('@ SAVE MODEL')
blast_model.save(out_path)
print('@ DONE')
    
