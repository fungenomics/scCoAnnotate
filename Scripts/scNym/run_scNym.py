#--------------- Libraries -------------------
import pandas as pd
from scnym.api import scnym_api
import anndata as ad
import sys
import scanpy as sc
import os
import random

# Set seed
random.seed(123456)

#def remove_batch(cell_id, batch):
#    return cell_id.replace(f'-{batch}', '')

def remove_batch(cell_id, batch):
    parts = cell_id.rsplit(f'-{batch}', 1)
    return ''.join(parts)

#--------------- Parameters -------------------
ref_path = str(sys.argv[1])
lab_path = str(sys.argv[2])
sample_path = str(sys.argv[3])
out_path = str(sys.argv[4])
threshold = float(sys.argv[5])
out_other_path = os.path.dirname(str(sys.argv[4]))

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

print('@ READ QUERY')
query = pd.read_csv(sample_path,
                    index_col=0,
                    sep=',',
                    engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

print('@ DONE')

#------------- Preparing input -------------

#Normalizing reference
ad_ref = ad.AnnData(X = ref,
                    obs = dict(obs_names=ref.index.astype(str),
                               label = labels['label']),
                    var = dict(var_names=ref.columns.astype(str))
                    )

# Now I normalize the matrix with scanpy:
# Normalize each cell by total counts over all genes,
# so that every cell has the same total count after normalization.
# If choosing `target_sum=1e6`, this is CPM normalization. 
sc.pp.normalize_total(ad_ref, target_sum=1e6)

# Logarithmize the data:
sc.pp.log1p(ad_ref)
sc.pp.filter_genes(ad_ref, min_cells=10)

# Normalizing query
ad_query = ad.AnnData(X = query,
                      obs = dict(obs_names=query.index.astype(str)),
                      var = dict(var_names=query.columns.astype(str))
)

# Now I normalize the matrix with scanpy:
# Normalize each cell by total counts over all genes,
# so that every cell has the same total count after normalization.
# If choosing `target_sum=1e6`, this is CPM normalization
sc.pp.normalize_total(ad_query, target_sum=1e6)

# Logarithmize the data:
sc.pp.log1p(ad_query)
sc.pp.filter_genes(ad_query, min_cells=10)

# Any cell with the annotation "Unlabeled" will be treated as part of the target 
# dataset and used for semi-supervised and adversarial training.
ad_query.obs["label"] = "Unlabeled"


#------------- Train scNym -------------

# Contatenate
adata = ad_ref.concatenate(ad_query)
file_model = out_other_path + '/model'
scnym_api(
    adata=adata,
    task='train',
    groupby='label',
    out_path=file_model,
    config='new_identity_discovery',
)

file_pred = out_other_path + '/predict'
scnym_api(
    adata=adata,
    task='predict',
    key_added='scNym',
    trained_model=file_model,
    out_path=file_pred,
    config='new_identity_discovery',
)

adata.obs['pred_label_reject'] = adata.obs.apply(lambda row: 'Unknown' if row['scNym_confidence'] < threshold else row['scNym'], axis=1)


adata.obs['cell_id'] = adata.obs.index.map(lambda cell_id: remove_batch(cell_id, adata.obs.loc[cell_id, 'batch']))

# Get only the query 
df = adata.obs[adata.obs["label"] == "Unlabeled"]
print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'cell': df.cell_id, "scNym": df.pred_label_reject})
pred_df.to_csv(out_path, index = False)
print('@ DONE')

#------------- Other outputs --------------

# Save data.frame output
print('@ WRITTING DATA FRAME OUTPUT ')
filename = out_other_path + '/whole_df_output.csv'
df.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')

# Save map object that contains the cell x label prob
print('@ WRITTING OUTPUT OBJECT ')
filename = out_other_path + '/scNym_output_object.h5ad'
adata.write_h5ad(filename= filename,
                  compression='gzip')
print('@ DONE ')

# make prob output matrix
pred_df['prob'] = df.scNym_confidence
pred_df = pred_df.pivot_table(index=pred_df.columns[0], columns='scNym', values='prob', fill_value=0).reset_index()
    
# rename column names 
pred_df.columns.name = None  # Remove the columns' name to match the R code
pred_df.columns = [''] + list(pred_df.columns[1:])

# save binary matrix
pred_df.to_csv(out_other_path + '/scNym_pred_score.csv', index=False)
