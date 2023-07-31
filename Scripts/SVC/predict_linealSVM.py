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
### Set seed
random.seed(123456) 

#--------------- Parameters -------------------
sample_path = str(sys.args[1])
model_path = str(sys.args[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.args[3]))
threshold = float(sys.args[4])
threads = int(sys.argv[5])
#--------------- Data -------------------------
print('@ READ QUERY')
query = pd.read_csv(sample_path,
                    index_col=0,
                    sep=',',
                    engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

print('@ DONE')

# load model 
print('@ LOAD MODEL')
SVM_model = pickle.load(open(model_path, 'rb'))
print('@ DONE')
### If the trained model was generated with a diff number of threads and the
## prediction is done with other number
SVM_model.n_jobs = threads
## We don't need to calculate again the predict_proba
#pred = SVM_model.predict(query.X)
pred_proba = SVM_model.predict_proba(query.X)
df = pd.DataFrame(pred_proba,index = query.obs_names,columns = SVM_model.classes_)
# Function to get the column name and maximum value for each row
def get_max_column_and_value(row):
    pred_label = row.idxmax()
    proba_label = row.max()
    return pred_label, proba_label

# Create a new column 'max_column' with the column name containing the maximum value for each row
df['pred_label'], df['proba_label'] = zip(*df.apply(get_max_column_and_value, axis=1))

#mm1 = pred.tolist()
#mm = df.max_column.tolist()
# mm1 == mm True
## It's the same

# Create a new column 'unknown_max_column' to store 'max_column' as 'unknown' if 'max_value' is lower than the threshold
df['pred_label_reject'] = df.apply(lambda row: 'Unknown' if row['proba_label'] < threshold else row['pred_label'], axis=1)

print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'cell': df.index, 'SVM': df.pred_label_reject})
pred_df.to_csv(out_path)
print('@ DONE')

#------------- Other outputs --------------
### Save the prob matrix
print('@ WRITTING PROB MATRIX ')
filename = out_other_path + '/prob_matrix.csv'
df.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')
