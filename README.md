# scCoAnnotate <img src ="https://user-images.githubusercontent.com/59002771/130340419-3d1eff0b-ecb2-4104-9bf4-1bb968aff433.png" width="50" height="50">

# Summary

Reference based prediction of cell-types using a fast and efficient Snakemake pipeline to increase automation and reduce the need to run several scripts and experiments. The pipeline allows the user to select what single-cell annotation tools they want to run on a selected reference to annotate a list of query datasets. It then outputs a consensus of the predictions across tools selected. This pipeline trains classifiers on genes common to the reference and all query datasets. 


# Installation and Dependencies

Tested with [R](https://www.r-project.org/) Version 4.2.2 and Python 3.11.2.

**R packages CRAN**
- Seurat
- tidyverse
- MetBrewer
- plotly
- caret
- Matrix

```
install.package(<pkg>)
```

Older version of Matrix package had to be installed for Seurat to work: https://github.com/satijalab/seurat/issues/6746

```
devtools::install_version("Matrix", version = "1.5.3", repos = "http://cran.us.r-project.org")
```

**R packages bioconductor**
- SingleCellExperiment
- SummarizedExperiment
- ComplexHeatmap
- WGCNA
- SingleR
- scClassify
- scuttle
- scran
- M3Drop

```
BiocManager::install(<pkg>)
```

**R packages github**
- singleCellNet (pcahan1/singleCellNet)
- scPred (powellgenomicslab/scPred)
- scibet (PaulingLiu/scibet)
- scLearn (bm2-lab/scLearn)
  
```
devtools::install_github(<pkg>)
```

**Python modules**
- numpy
- pandas
- scHPL
- sklearn
- anndata
- matplotlib
- scanpy
- datetime
- tensorflow
- tables
- snakemake

```
pip install <pkg>
```

# Quickstart

UPDATE 

##  Config File:

```yaml 
# UPDATE 
```

### An Example Config is attached 

```yaml 
# UPDATE
```

## Submission File:

An example of the submission file is also available in this repository and is called submit.sh. This is for TORQUE schedulers.


``` bash 
# UPDATE
```

# Tools Available
      - scPred
      - SingleR
      - scClassify
      - SciBet
      - singleCellNet
      - scHPL 
      - SVMlinear
      - scLearn
      - Correlation

# Adding New Tools:

``` R
# UPDATE 
```

# Tools 

## scClassify

Detailed documentation for scClassify train and predict scripts, written July 2023 by Bhavyaa Chandarana

* scCoAnnotate input reference and query have cells as the rows, genes as columns. scClassify (and the Seurat function used for normalization, see below) requires genes on the rows and cells on the columns. Therefore, I used `WGCNA::transposeBigData()` (function optimized for large sparse matrices) to transpose the inputs before normalization and training/prediction.

* scClassify documentation defines "log-transformed" data as "size-factor normalized" data ([source](https://www.bioconductor.org/packages/devel/bioc/vignettes/scClassify/inst/doc/scClassify.html#2_Setting_up_the_data)). Function documentation for both `train_scClassify()` and `predict_scClassify()` specify that reference and query data must be "log-transformed" ([source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf)) Therefore, I am normalizing both query and reference with `Seurat::NormalizeData()` (default parameters), which performs log normalization with scale factor 10000 ([source](https://satijalab.org/seurat/reference/normalizedata))

* scClassify train and predict functions `train_scClassify()` and `predict_scClassify()` both allow parallelization with package `BiocParallel`. If greater than one thread was requested by the user, I turn parallelization mode on with parallel = TRUE, and set the `BiocParallel` parameter to `BiocParallel::MulticoreParam()` with workers equal to number of requested threads (based on code in [this issue](https://github.com/SydneyBioX/scClassify/issues/14)) Otherwise I have set parallelization to FALSE and the `BiocParallel` parameter to `BiocParallel::SerialParam()` (which is the default value of the parameter in the functions - [source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf)).

* scClassify train function `train_scClassify()` can either return a list output for the model, or an R object of class `scClassifyTrainModel`, based on boolean argument `returnList`. Both types work as input for prediction with `predict_scClassify()`. However, in order to use `scClassify::cellTypeTree()` to extract and output the tree produced by scClassify during training, the input must be the R object of class `scClassifyTrainModel`. Therefore, I have chosen to set `returnList` in `train_scClassify()` to FALSE (default: TRUE), and use the resulting object for `cellTypeTree()`. (Function documentation [source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf))

* `scClassify::plotCellTypeTree()` produces a ggplot object. Therefore, I am using `ggplot2::ggsave()` to save it as a png file. (Function documentation [source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf))

## scPred

Both reference and query is normaluzed using `Seurat::NormalizeData()`.     
Needs computed PCA space. Dims set to 1:30 according to tutorial.    
Default model `SVMradial`. Option to switch model should be set up in snakemake.   

Normalization and parameters based on this tutorial:   
https://powellgenomicslab.github.io/scPred/articles/introduction.html

## SingleR 

Both reference and query is normaluzed using `scuttle::logNormCounts()`. Both reference and query is converted to SingleCellExperiment objects before normalization.   

Deviation from default parameters: 
* `de.method = de.method="wilcox"`
Method for generating marker genes for each class in reference. Wilcox is recomended when single cell data is used as reference

Normalization and parameters based on this tutorial:
http://www.bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html#3_Using_single-cell_references

## singleCellNet

Documentation written by: Rodrigo Lopez Gutierrez
Date written: 2023-08-01

Input for `singleCellNet` is raw counts for both reference and query. The reference is normalized within the `scn_train()` function. The query is currently not normalized. In the tutorial example they used raw query data. Furthermore, according to the tutorial, the classification step is robust to the normalization and transformation steps of the query data sets. They claim that one can even directly use raw data to query and still obtains accurate classification. This could be tested in the future with our data to see if normalized queries perform better.

Normal parameters were used in both the training and prediction functions, with the expection of the following parameters:
* In `scn_train()`, we used parameter `nTrees = 500` compared to the default `nTrees = 1000`. This parameter changes the number of trees for the random forest classifier. The value selected is based on Hussein's thesis and is changed to improve the speed of `singleCellNet`. It is mentioned that additional training parameters may need to be adjusted depending on the quality of the reference data. Additionally, tutorial mentions that classifier performance may increase if the values for `nTopGenes` and `nTopGenePairs` are increased.
* In `scn_predict()`, we used parameter `nrand = 0` compared to the default `nrand = 50`. This parameter changes the number of randomized single cell RNA-seq profiles which serve as positive controls that should be mostly classified as `rand` (unknown) category. If left at default value, then this would generate extra cells that might complicate downstream consolidation of the consensus predictions for each cell. Therefore, the selected value is used to avoid complication. 

singleCellNet workflow was generated following the tutorial below:
https://pcahan1.github.io/singleCellNet/

## Correlation

Documentation written by: Rodrigo Lopez Gutierrez
Date written: 2023-08-02

The Correlation tool runs a correlation-based cell type prediction on a sample of interest, given the mean gene expression per label for a reference.
The function to label by Spearman correlation was originally generated by Selin Jessa and Marie Coutlier
Path to original file: `/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/stable/code/scripts/predict_celltype_cor.R`

Input for `Correlation` is raw counts for both reference and query. Both the reference and the query are normalized using `Seurat::NormalizeData()`.

Training script generates a matrix with the mean gene expression for each label in the reference.
Prediction script calculates a correlation between each cell in the query and each label in mean gene expression matrix generated in the training script. Then we assign each cell the most highly correlated label. 
* `label_correlation()` function has a parameter `threshold_common_genes` which sets the percentage of query dataset genes required to be in the reference dataset in order to proceed. This parameter is currently not utilized as the preprocessing done in the beginning of the snakefile is extracting only the common genes between the reference and the queries.

Currently only outputting a table with each cell, the most highly correlated label, and the corresponding correlation score for that label. In the future we could export the full correlation matrix, if necessary.

## scLearn

Detailed documentation for scLearn train and predict scripts, written August 2023 by Bhavyaa Chandarana. Added information Tomas Vega Waichman in 2023-08-04.

Preprocessing performed in the same way as scLearn documentation tutorial [source](https://github.com/bm2-lab/scLearn#tutorial) for `Single-label single cell assignment`

* scCoAnnotate input reference and query have cells as the rows, genes as columns. scLearn requires genes on the rows and cells on the columns. Therefore, I used `WGCNA::transposeBigData()` (function optimized for large sparse matrices) to transpose the inputs before normalization and training/prediction.

* In order to avoid cell filtering, the reference and query matrix were normalized using Seurat::NormalizeData since the authors make the logNormalization in this way but manually (using a scale.factor = 10000 and then log(ref + 1)). Because of this the arguments of 'species' is not used and this allows to use this methods in order species different to the Human and Mouse.

* Used default value `10` for argument `bootstrap_times` in training function. According to tool documentation, this can be increased to improve accuracy for unassigned cells(?) but increase train time.

* Default parameters were used for tool prediction 

*  Added some outputs, for prediction added a table with the selected genes for the model. In prediction added and output with the whole data.frame with the probabilities for each cell.

## singleCellNet

Documentation written by: Rodrigo Lopez Gutierrez
Date written: 2023-08-01

Input for `singleCellNet` is raw counts for both reference and query. The reference is normalized within the `scn_train()` function. The query is currently not normalized. In the tutorial example they used raw query data. Furthermore, according to the tutorial, the classification step is robust to the normalization and transformation steps of the query data sets. They claim that one can even directly use raw data to query and still obtains accurate classification. This could be tested in the future with our data to see if normalized queries perform better.

Normal parameters were used in both the training and prediction functions, with the expection of the following parameters:
* In `scn_train()`, we used parameter `nTrees = 500` compared to the default `nTrees = 1000`. This parameter changes the number of trees for the random forest classifier. The value selected is based on Hussein's thesis and is changed to improve the speed of `singleCellNet`. It is mentioned that additional training parameters may need to be adjusted depending on the quality of the reference data. Additionally, tutorial mentions that classifier performance may increase if the values for `nTopGenes` and `nTopGenePairs` are increased.
* In `scn_predict()`, we used parameter `nrand = 0` compared to the default `nrand = 50`. This parameter changes the number of randomized single cell RNA-seq profiles which serve as positive controls that should be mostly classified as `rand` (unknown) category. If left at default value, then this would generate extra cells that might complicate downstream consolidation of the consensus predictions for each cell. Therefore, the selected value is used to avoid complication. 

singleCellNet workflow was generated following the tutorial below:
https://pcahan1.github.io/singleCellNet/
