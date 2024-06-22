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
sample_path = str(sys.argv[1])
model_path = os.path.dirname(str(sys.argv[2]))
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threshold = float(sys.argv[4])
threads = int(sys.argv[5])

cb.config.N_JOBS = threads
cb.config.RANDOM_SEED = 123456

#--------------- Data -------------------

# read query matrix
print('@ READ QUERY')
query = pd.read_csv(sample_path,
                    index_col=0,
                    sep=',',
                    engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

print('@ DONE')

# Query preprocessing
query = ad.AnnData(X = query,
                   obs = dict(obs_names=query.index.astype(str)),
                   var = dict(var_names=query.columns.astype(str))
                   )

print('@ LOADING MODEL')
blast_model = cb.blast.BLAST.load(model_path)
print('@ DONE')

#----------- Predict CellBlast --------
query_hits = blast_model.query(query)

query_hits = query_hits.reconcile_models().filter(by="pval", cutoff=0.05)

predictions = query_hits.annotate(field = "labels",
                                  majority_threshold = threshold)

print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'cell': predictions.index,
                       'CellBlast': predictions.labels})

pred_df.to_csv(out_path,
               index = False)
print('@ DONE')

#----------------------------------------
### The prob only return the probabilities of the label assigned, so it's not usefull for CAWPE.
## Then we keep doing the same procedure of making the binary matrix.
# make binary output matrix
pred_df['prob'] = 1
pred_df = pred_df.pivot_table(index=pred_df.columns[0], columns='CellBlast', values='prob', fill_value=0).reset_index()
    
# rename column names 
pred_df.columns.name = None  # Remove the columns' name to match the R code
pred_df.columns = [''] + list(pred_df.columns[1:])

# save binary matrix
pred_df.to_csv(out_other_path + '/CellBlast_pred_score.csv', index=False)
