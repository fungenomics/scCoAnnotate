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
from scHPL import learn, train, utils
from scHPL.utils import TreeNode, _print_node, _count_nodes
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
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
classifier = str(sys.argv[4])
dimred = bool(sys.argv[5])

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

# A tree could be use an input, in that case the read_tree can be used.
# The tree format is Newick formatted file
# tree = utils.read_tree(treePath)

# They provided this solution here: https://github.com/lcmmichielsen/scHPL/issues/7
tree = utils.create_tree('root')
tree = learn._construct_tree(tree, labels)

#------------- Train scHPL -------------
# classifier could be svm, svm_occ, knn
# - the linear SVM when your integrated data still has a lot of
# dimensions (e.g. when you have used Seurat to integrate the datasets)
# - the kNN when your integrated data has less, 10-50, dimensions
# (e.g. when you have used scVI or Harmony to integrate the datasets)
# - the one-class SVM when your main focus is to find unseen cell
# populations. A downside of the one-class SVM, however, is that the
# classification performance drops.
tree = train.train_tree(adata.X,
                        labels.label,
                        tree,
                        classifier = classifier,
                        dimred = dimred, 
                        useRE = True,
                        FN = 0.5)

# Save the model to disk
print('@ SAVE MODEL')
pickle.dump(tree, open(out_path, 'wb'))
print('@ DONE')


#------------- Other outputs --------------
# Plot the tree
# I'm using this method since they provided it here:
# https://github.com/lcmmichielsen/scHPL/issues/5
def _print_node(node, hor, ver_steps, fig, new_nodes):
    global ver
    # Add horizontal line
    x, y = ([np.max([0.05, hor-0.045]), hor], [ver, ver])
    line = mlines.Line2D(x,y, lw=1)
    fig.add_artist(line)
    
    # Add textbox
    if np.isin(node.name[0], new_nodes):
        txt = r"$\bf{" + node.name[0] + "}$"
    else:
        txt = node.name[0]
    
    for n in node.name:
        if(n != node.name[0]):
            if np.isin(n, new_nodes):
                txt = txt + ' & ' + r"$\bf{" + n + "}$"
            else:
                txt = txt + ' & ' + n
                
    fig.text(hor,ver, txt, size=10,
             ha = 'left', va='center',
             bbox = dict(boxstyle='round', fc='w', ec='k'))
    
    # Continue with child nodes
    hor = hor+0.05
    ver_line_start = ver
    ver_line_end = ver
    
    for i in node.descendants:
        ver = ver-ver_steps
        ver_line_end = ver
        _print_node(i, hor, ver_steps, fig, new_nodes)
        
    # Add vertical line
    x, y = ([np.max([0.05, hor-0.045]), np.max([0.05, hor-0.045])], 
            [ver_line_start, ver_line_end])
    line = mlines.Line2D(x,y, lw=1)
    fig.add_artist(line)
    
def print_tree(filepath,
               tree: TreeNode, 
               new_nodes: list = []):
    '''Print the tree
        Parameters
        ----------
        tree : TreeNode
            Tree to print
        new_nodes : List = []
            Nodes recently added to the tree, these are printed in bold
    
        Returns
        -------
        None.
    '''
    
    global ver
    ver = 0.93
    
    count = _count_nodes(tree)
    ver_steps = 0.9/count
    plot_height = count*0.3
    fig = plt.figure(figsize=(6,plot_height)) # This size is hard coded
    ax = plt.subplot(111)

    _print_node(tree[0], hor=0.05, ver_steps=ver_steps, fig=fig, 
                new_nodes = new_nodes)
    
    plt.axis('off')
    plt.savefig(filepath, dpi=1000)
    plt.close()
    
print('@ PLOTTING')
filepath = out_other_path + '/tree.png'
print_tree(filepath=filepath,
           tree = tree)
print('@ DONE')
