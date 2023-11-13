#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 17:52:32 2023

@author: tomas.vega
"""
## Python 3.11.2
#--------------- Libraries -------------------
import pandas as pd
from scHPL import predict
import anndata as ad
import sys
import os
import scanpy as sc
import pickle
import random

# Set seed
random.seed(123456) 

#--------------- Parameters -------------------
sample_path = str(sys.argv[1])
model_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threshold = float(sys.argv[4])
#threads = as.numeric(args[4])

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
scHPL_model = pickle.load(open(model_path, 'rb'))
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

#----------- Predict scHPL --------

pred = predict.predict_labels(testdata= query.X,
                              tree = scHPL_model,
                              threshold = threshold)

print('@ WRITTING PREDICTIONS')

pred_labels = pd.DataFrame({'cell': query.obs_names, 'scHPL': pred})
pred_labels.to_csv(out_path, index = False)
print('@ DONE')
#----------------------------------------
