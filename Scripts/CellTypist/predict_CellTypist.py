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
sample_path = str(sys.argv[1])
model_path = str(sys.argv[2])
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threshold = float(sys.argv[4])
threads = int(sys.argv[5])
majority_voting = bool(sys.argv[6])

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

## Now I normalize the matrix with scanpy:
#Normalize each cell by total counts over all genes,
#so that every cell has the same total count after normalization.
#If choosing `target_sum=1e6`, this is CPM normalization
#1e4 similar as Seurat
sc.pp.normalize_total(query, target_sum=1e4)
#Logarithmize the data:
sc.pp.log1p(query)


# load model 
print('@ LOAD MODEL'))
CellTypist_model = pickle.load(open(model_path, 'rb'))
print('@ DONE')

### If the trained model was generated with a diff number of threads and the
## prediction is done with other number
CellTypist_model.classifier.n_jobs = threads

#------------- Predict CellTypist -------------

predictions = celltypist.annotate(query,
                                  model = CellTypist_model,
                                  majority_voting = majority_voting,
                                  mode = 'prob match',
                                  p_thres = threshold)

## If the majority voting is true, return that result
if majority_voting:
  pred_df = pd.DataFrame({'cell': predictions.predicted_labels.index,
                          'CellTypist': predictions.predicted_labels.majority_voting})
else:
  pred_df = pd.DataFrame({'cell': predictions.predicted_labels.index,
                        "CellTypist": predictions.predicted_labels.predicted_labels})

print('@ WRITTING PREDICTIONS')
pred_df.to_csv(out_path, index = False)
print('@ DONE')

#------------- Other outputs --------------
## Save the outputs as .csv
print('@ WRITTING CSV OUTPUTS')
predictions.to_table(out_other_path)
print('@ DONE')
## Save the output plots
print('@ GENERATING OUTPUT PLOTS')
prediction.to_plot(out_other_path)
print('@ DONE')

### Save the prob matrix
print('@ WRITTING OUTPUT OBJECT')
filename = out_other_path + '/CellTypist_output_object.h5ad'
adata = predictions.to_adata(insert_prob = True)
adata.write_h5ad(filename= filename,
                  compression='gzip')
print('@ DONE ')

                                  
