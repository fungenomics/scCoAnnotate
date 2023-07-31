# scCoAnnotate <img src ="https://user-images.githubusercontent.com/59002771/130340419-3d1eff0b-ecb2-4104-9bf4-1bb968aff433.png" width="50" height="50">

scRNA-seq based prediction of cell-types using an automated Snakemake pipeline to produce a consensus of several prediction tools. 

# Summary

The pipeline allows the user to select what single-cell annotation tools they want to run on a selected reference to annotate a list of query datasets. It then outputs a consensus of the predictions across tools selected. This pipeline trains classifiers on genes common to the reference and all query datasets. The pipeline also features parallelization options to exploit the computational resources available. 

# Installation and Dependencies

* Install or load [R](https://www.r-project.org/) version 4.2.2 and Python version 3.6.5.

* Install [Snakemake](https://snakemake.readthedocs.io/en/stable/) version 5.32.0 in your linux environment.

```bash
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

* Install all the dependencies for the tools you plan on using. See [the list of dependencies below](#Required-packages)
  
* Clone this repository into your directory of choice

# Quickstart

The processes of the snakemake pipeline are arranged as per this rule graph:

![dag](https://user-images.githubusercontent.com/59002771/191146873-5c680bbd-d11c-418c-ae96-7662ee7f99ed.png)

Rule preprocess gets the common genes and creates temporary reference and query datasets based ob the common genes. Rule concat appends all predictions into one tab seperate file (`prediction_summary.tsv`) and gets the consensus prediction.

After preparing your config file, you can run the pipeline with the following command:

```bash
snakemake --use-conda --configfile config.yml --cores 3
```

##  Config file template
```yaml 
# target directory
output_dir: <path to outputs directory>
# path to reference to train classifiers on (cell x gene raw counts)
training_reference: <path to reference csv file with RAW counts per cell, genes as columns and cells as rows>
# path to annotations for the reference (csv file with cellname and label headers)
reference_annotations: <csv with labels per cell, the column header for the labels should be "label">
# path to query datasets (cell x gene raw counts)
query_datasets:
      - <path to query csv file 1 with RAW counts per cell, genes as columns and cells as rows>
      - <path to query csv file 2 with RAW counts per cell, genes as columns and cells as rows>
      .
      .
# step to check if required genes are kept between query and reference
check_genes: False
# path for the genes required
genes_required: Null
# rejection option for SVM
rejection: True
# classifiers to run
tools_to_run:
      - <tool 1>
      - <tool 2>
      .
      .
# benchmark tools on reference
benchmark: False
plots: True
consensus:
      - <tool 1>
      - <tool 2>
      .

```


### Example config file 

```yaml 
# target directory
output_dir: /project/kleinman/hussein.lakkis/from_hydra/test
# path to reference to train classifiers on (cell x gene raw counts)
training_reference: /project/kleinman/hussein.lakkis/from_hydra/2022_01_10-Datasets/Cortex/p0/expression.csv
# path to annotations for the reference (csv file with cellname and label headers)
reference_annotations: /project/kleinman/hussein.lakkis/from_hydra/2022_01_10-Datasets/Cortex/p0/labels.csv
# path to query datasets (cell x gene raw counts)
query_datasets:
      - /project/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/BT2016062/expression.csv
      - /project/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/P-1694_S-1694_multiome/expression.csv
      - /project/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/P-1701_S-1701_multiome/expression.csv
# step to check if required genes are kept between query and reference
check_genes: False
# path for the genes required
genes_required: Null
# rejection option for SVM
rejection: True
# classifiers to run
tools_to_run:
      - ACTINN
      - scHPL
      - scClassify
      - correlation
      - scmapcluster
      - scPred
      - SingleCellNet
      - SVM_reject
      - SingleR
      - CHETAH
      - scmapcell
      - SciBet
# benchmark tools on reference
benchmark: False
plots: True
consensus:
      - all
```

## Submission File:

An example of the submission file is also available in this repository and is called submit.sh. This is for TORQUE schedulers.

``` bash 
#!/usr/bin/bash
#PBS -N scCoAnnotate
#PBS -o logs/err.txt
#PBS -e logs/out.txt
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=3
#PBS -l mem=125G
#PBS -l vmem=125G

# you need to do this step 

cd {root directory of the pipeline}
mkdir -p logs

# set up the environment
#conda init bash
module load miniconda/3.8
source ~/.conda_init.sh
module load snakemake/5.32.0
module load hydra/R/4.0.5
module load python/3.6.5

# Run the snakemake
snakemake --use-conda --configfile config.yml --cores 3
```

# Available tools

1. [ACTINN](https://github.com/mafeiyang/ACTINN)
2. [SciBet](https://github.com/PaulingLiu/scibet)
4. [Spearman Correlation](https://statistics.laerd.com/statistical-guides/spearmans-rank-order-correlation-statistical-guide.php)
5. [SVM](https://scikit-learn.org/stable/modules/svm.html)
6. SVM Rejection
7. [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
8. [SingleCellNet](https://github.com/pcahan1/singleCellNet)
9. [CHETAH](https://www.bioconductor.org/packages/release/bioc/html/CHETAH.html)
10. [scHPL](https://github.com/lcmmichielsen/scHPL)
11. [scPred](https://github.com/powellgenomicslab/scPred)
12. [scmap (cell and cluster)](https://bioconductor.org/packages/release/bioc/html/scmap.html)

# Required packages

## Python Libraries:

```
tensorboard==1.7.0
tensorboard-data-server==0.6.1
tensorboard-plugin-wit==1.8.0
tensorflow==1.7.0
tensorflow-estimator==2.5.0
sklearn==0.0
scikit-learn==0.24.1
pandas==1.1.5
numpy==1.19.5
numpy-groupies==0.9.13
numpydoc==1.1.0
scHPL==0.0.2
```

## R Libraries:

```
glue
data.table
scPred_1.9.2
SingleCellExperiment_1.12.0
SummarizedExperiment_1.20.0
CHETAH_1.6.0
scmap_1.12.0 
singleCellNet == 0.1.0
scibet == 1.0
SingleR == 1.4.1
Seurat == 4.0.3
dplyr == 1.0.7
tidyr == 1.1.3
viridis == 0.6.1
ggsci == 2.9
tidyverse == 1.3.1
```
# Adding New Tools:

to add new tools, you have to add this template to the the snakefile as such:

``` R
rule {rulename}:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/{rulename}/{rulename}_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/{rulename}/{rulename}_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/{rulename}/{rulename}_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/{rulename}/{rulename}.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_{rulename}.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"   

 ```   
 The tool script you add must generate outputs that match the output of the rule.

# scClassify

Detailed documentation for scClassify train and predict scripts, written July 2023 by Bhavyaa Chandarana

* scCoAnnotate input reference and query have cells as the rows, genes as columns. scClassify (and the Seurat function used for normalization, see below) requires genes on the rows and cells on the columns. Therefore, I used `WGCNA::transposeBigData()` (function optimized for large sparse matrices) to transpose the inputs before normalization and training/prediction.

* scClassify documentation defines "log-transformed" data as "size-factor normalized" data ([source](https://www.bioconductor.org/packages/devel/bioc/vignettes/scClassify/inst/doc/scClassify.html#2_Setting_up_the_data)). Function documentation for both `train_scClassify()` and `predict_scClassify()` specify that reference and query data must be "log-transformed" ([source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf)) Therefore, I am normalizing both query and reference with `Seurat::NormalizeData()` (default parameters), which performs log normalization with scale factor 10000 ([source](https://satijalab.org/seurat/reference/normalizedata))

* scClassify train and predict functions `train_scClassify()` and `predict_scClassify()` both allow parallelization with package `BiocParallel`. If greater than one thread was requested by the user, I turn parallelization mode on with parallel = TRUE, and set the `BiocParallel` parameter to `BiocParallel::MulticoreParam()` with workers equal to number of requested threads (based on code in [this issue](https://github.com/SydneyBioX/scClassify/issues/14)) Otherwise I have set parallelization to FALSE and the `BiocParallel` parameter to `BiocParallel::SerialParam()` (which is the default value of the parameter in the functions - [source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf)).

* scClassify train function `train_scClassify()` can either return a list output for the model, or an R object of class `scClassifyTrainModel`, based on boolean argument `returnList`. Both types work as input for prediction with `predict_scClassify()`. However, in order to use `scClassify::cellTypeTree()` to extract and output the tree produced by scClassify during training, the input must be the R object of class `scClassifyTrainModel`. Therefore, I have chosen to set `returnList` in `train_scClassify()` to FALSE (default: TRUE), and use the resulting object for `cellTypeTree()`. (Function documentation [source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf))

* `scClassify::plotCellTypeTree()` produces a ggplot object. Therefore, I am using `ggplot2::ggsave()` to save it as a png file. (Function documentation [source](https://www.bioconductor.org/packages/release/bioc/manuals/scClassify/man/scClassify.pdf))

# scPred

Both reference and query is normaluzed using `Seurat::NormalizeData()`.     
Needs computed PCA space. Dims set to 1:30 according to tutorial.    
Default model `SVMradial`. Option to switch model should be set up in snakemake.   

Normalization and parameters based on this tutorial:   
https://powellgenomicslab.github.io/scPred/articles/introduction.html

# SingleR 

Both reference and query is normaluzed using `scuttle::logNormCounts()`. Both reference and query is converted to SingleCellExperiment objects before normalization.   

Deviation from default parameters: 
* `de.method = de.method="wilcox"`
Method for generating marker genes for each class in reference. Wilcox is recomended when single cell data is used as reference

Normalization and parameters based on this tutorial:    
http://www.bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html#3_Using_single-cell_references
