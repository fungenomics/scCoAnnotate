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
from sklearn.svm import SVC
from sklearn.calibration import CalibratedClassifierCV
import anndata as ad
import sys
import scanpy as sc
import pickle
import os
import random

# Set seed
random.seed(123456) 

# Function to get the column name and maximum value for each row
def get_max_column_and_value(row):
    pred_label = row.idxmax()
    proba_label = row.max()
    return pred_label, proba_label

#--------------- Parameters -------------------

sample_path = str(sys.argv[1])
model_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threshold = float(sys.argv[4])
threads = int(sys.argv[5])
tool_name = str(sys.argv[6])

#--------------- Data -------------------------

print('@ READ QUERY')
query = pd.read_csv(sample_path,
                    index_col=0,
                    sep=',',
                    engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

print('@ DONE')

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

#Logarithmize the data:
sc.pp.log1p(query)

# load model 
print('@ LOAD MODEL')
SVM_model = pickle.load(open(model_path, 'rb'))
print('@ DONE')


# If the trained model was generated with a diff number of threads and the
# prediction is done with other number
SVM_model.n_jobs = threads

# Calculate predicted probability
pred_proba = SVM_model.predict_proba(query.X)
df = pd.DataFrame(pred_proba,index = query.obs_names,columns = SVM_model.classes_)

# Create a new column 'max_column' with the column name containing the maximum value for each row
df['pred_label'], df['proba_label'] = zip(*df.apply(get_max_column_and_value, axis=1))

# Create a new column 'unknown_max_column' to store 'max_column' as 'unknown' if 'max_value' is lower than the threshold
df['pred_label_reject'] = df.apply(lambda row: 'Unknown' if row['proba_label'] < threshold else row['pred_label'], axis=1)

print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'cell': df.index, tool_name: df.pred_label_reject})
pred_df.to_csv(out_path, index = False)
print('@ DONE')

#------------- Other outputs --------------

# Save the prob matrix
print('@ WRITTING PROB MATRIX ')
filename = out_other_path + '/prob_matrix.csv'
df.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')
