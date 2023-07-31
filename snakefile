# import libraries
import os

# Get the names of query samples from the paths given in the query section of the config
samples = [os.path.basename(os.path.dirname(query_path)) for query_path in config['query_datasets']]

# Create subdirectories for each query sample 
for sample in samples:
  path = "{output_dir}/{sample}".format(output_dir = config["output_dir"], sample = sample)
  if not os.path.exists(path):
    os.mkdir(path)

"""
One rule that directs the DAG which is represented in the rulegraph
"""
rule all:
  input:
    expand(config["output_dir"] + "/{sample}/Prediction_Summary.tsv", 
           sample  = samples),
    expand(config["output_dir"] + "/{sample}/{tool}/{tool}_pred.csv", 
           tool = config['tools_to_run'],
           sample = samples)

"""
rule that gets the gets the interesction in genes between samples and reference
It outputs temporary reference and query datasets based on the common genes
"""
rule preprocess:
  input:
    reference = config['training_reference'],
    query = config['query_datasets']
  output:
    reference = config['output_dir'] + "/expression.csv",
    query = expand(config['output_dir'] + "/{sample}/expression.csv", sample = samples)
  params:
    check_genes = bool(config['check_genes']),
    genes_required = config['genes_required']
  log: config['output_dir'] + "/preprocess.log"
  priority: 50
  shell:
    "Rscript {workflow.basedir}/Scripts/preprocess.R "
    "--ref {input.reference} "
    "--query {input.query} "
    "--output_dir {config[output_dir]} "
    "--check_genes {params.check_genes} "
    "--genes_required {params.genes_required} "
    "&> {log} "

"""
rule that gets the Consensus 
"""
rule concat:
  input:
    results = expand(config["output_dir"] + "/{sample}/{tool}/{tool}_pred.csv",
                     tool=config['tools_to_run'],
                     sample  = samples),
    sample = expand(config["output_dir"] + "/{sample}",
                    sample  = samples)
  output:
    expand(config["output_dir"] + "/{sample}/Prediction_Summary.tsv", 
           sample = samples)
  log: 
    expand(config["output_dir"] + "/{sample}/Gatherpreds.log", sample = samples)
  shell:
    "python3 {workflow.basedir}/Scripts/Gather_Preds.py "
    "-i {input.sample} "
    "-c {config[consensus]} "
    "-k {input.sample} "
    "&> {log}"    

"""
Rules for R based tools.


#---------------------------------------------------------------------------

rule train_SingleR:
  input:
    reference = config['output_dir'] + "/expression.csv",
    labfile = config['reference_annotations']
  output:
    model = config['output_dir'] + "/SingleR/SingleR_model.Rda"
  params:
      basedir = {workflow.basedir}
  log: 
    config['output_dir'] + "/SingleR/SingleR.log"
  benchmark:
    config['output_dir'] + "/SingleR/SingleR_train_benchmark.txt"
  threads: 1
  resources: 
  shell:
    """
    Rscript {params.basedir}/Scripts/SingleR/train_SingleR.R \
    {input.reference} \
    {input.labfile} \
    {output.model} \
    {threads} \
    &> {log}
    """

rule predict_SingleR:
  input:
    query = config['output_dir'] + "/{sample}/expression.csv",
    model = config['output_dir'] + "/SingleR/SingleR_model.Rda"
  output:
    pred = config['output_dir'] + "/{sample}/SingleR/SingleR_pred.csv"
  params:
      basedir = {workflow.basedir}
  log: 
    config['output_dir'] + "/{sample}/SingleR/SingleR.log"
  benchmark:
    config['output_dir'] + "/{sample}/SingleR/SingleR_predict_benchmark.txt"
  threads: 1
  resources: 
  shell:
    """
    Rscript {params.basedir}/Scripts/SingleR/predict_SingleR.R \
    {input.query} \
    {input.model} \
    {output.pred} \
    {threads} \
    &> {log}
    """

rule train_scPred:
  input:
    reference = config['output_dir'] + "/expression.csv",
    labfile = config['reference_annotations']
  output:
    model_type = config['output_dir'] + "/scPred/scPred_model.Rda"
  params:
    basedir = {workflow.basedir},
    model = "svmRadial"
  log: 
    config['output_dir'] + "/scPred/scPred.log"
  benchmark:
    config['output_dir'] + "/scPred/scPred_train_benchmark.txt"
  threads: 1
  resources: 
  shell:
    """
    Rscript {params.basedir}/Scripts/scPred/train_scPred.R \
    {input.reference} \
    {input.labfile} \
    {output.model} \
    {threads} \
    {params.model_type} \
    &> {log}
    """

rule predict_scPred:
  input:
    query = config['output_dir'] + "/{sample}/expression.csv",
    model = config['output_dir'] + "/scPred/scPred_model.Rda"
  output:
    pred = config['output_dir'] + "/{sample}/scPred/scPred_pred.csv"
  params:
      basedir = {workflow.basedir}
  log: 
    config['output_dir'] + "/{sample}/scPred/scPred.log"
  benchmark:
    config['output_dir'] + "/{sample}/scPred/scPred_predict_benchmark.txt"
  threads: 1
  resources: 
  shell:
    """
    Rscript {params.basedir}/Scripts/scPred/predict_scPred.R \
    {input.query} \
    {input.model} \
    {output.pred} \
    {threads} \
    &> {log}
    """

#---------------------------------------------------------------------------

rule correlation:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/correlation/correlation_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/correlation/correlation_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/correlation/correlation_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/correlation/correlation.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript {workflow.basedir}/Scripts/run_correlation.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule SciBet:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/SciBet/SciBet_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SciBet/SciBet_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SciBet/SciBet_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SciBet/SciBet.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SciBet.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"
 
rule SingleCellNet:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SingleCellNet.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"   
    
rule CHETAH:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/CHETAH/CHETAH_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    total_time = expand("{output_dir}/{sample}/CHETAH/CHETAH_total_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/CHETAH/CHETAH.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_CHETAH.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"
    
rule scmapcluster:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/scmapcluster/scmapcluster.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scmapcluster.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule scmapcell:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scmapcell/scmapcell_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scmapcell/scmapcell_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scmapcell/scmapcell_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/scmapcell/scmapcell.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scmapcell.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

"""
Rules for python based tools.
"""

rule ACTINN:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/ACTINN/ACTINN_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/ACTINN/ACTINN_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/ACTINN/ACTINN_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/ACTINN/ACTINN.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_ACTINN.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}" 
    
rule SVM_reject:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  params:
    rejection = config.get("rejection", "")
  output:
    pred = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SVM_reject/SVM_reject.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_SVM_reject.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
        "--rej {params.rejection}   "
    "--output_dir {input.output_dir} "
    "&> {log}"  

rule scHPL:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/scHPL/scHPL_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scHPL/scHPL_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scHPL/scHPL_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scHPL/scHPL.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_scHPL.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"  

rule scPred:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scPred/scPred_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scPred/scPred_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scPred/scPred_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scPred/scPred.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scPred.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"   

rule SingleR:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/SingleR/SingleR_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SingleR/SingleR_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SingleR/SingleR_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/SingleR/SingleR.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SingleR.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"
