#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 18:37:59 2023

@author: tomas.vega
"""
## Python 3.11.2

#--------------- Libraries -------------------
import numpy as np
import pandas as pd
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import anndata as ad
import sys
import scanpy as sc
import pickle
import os
import random

# Set seed
random.seed(123456) 

#--------------- Parameters -------------------
ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threads = int(sys.argv[4])

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
                   obs = dict(obs_names=ref.index.astype(str)),
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

## Transform labels to array

label = np.array(labels['label'])

#------------- Train SVM Lineal -------------
# kernel could be ‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’
# When the constructor option probability is set to True, class membershipsq
# probability estimates (from the methods predict_proba and predict_log_proba)
# are enabled. 
svm_model = LinearSVC()

# Calibrate the model using 5-fold cross validation
calibrated_model = CalibratedClassifierCV(svm_model,
                                          n_jobs = threads) #Default
calibrated_model.fit(adata.X, label)

#Save the model to disk
print('@ SAVE MODEL')
pickle.dump(calibrated_model, open(out_path, 'wb'))
print('@ DONE')

#------------- Other outputs --------------
# Plot the tree
filepath = out_other_path + '/train_parameters.csv'
print('@ WRITE TRAIN PARAMETERS ')
params_table = calibrated_model.get_params()
params_table = pd.DataFrame.from_dict(params_table,
                                      orient= 'index')
params_table.to_csv(filepath)
print('@ DONE')

# I cannot get the features per class, since:
# there is attribute coef_ for SVM classifier but it only works for SVM with
# linear kernel. For other kernels it is not possible because data are transformed
# by kernel method to another space, which is not related to input space
# (https://stackoverflow.com/questions/21260691/how-to-obtain-features-weights)
# I cannot get the scores since I use the CalibratedClassifierCV too, 
# and it doens't have the score attribute.
