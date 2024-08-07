
# set parameters related to converting gene names 
def set_gene_conversion_parameters(config):
  for ref in config['references'].keys():
    try:
      config["references"][ref]['convert_ref_mm_to_hg']
    except:
      config["references"][ref]['convert_ref_mm_to_hg'] = False

# set parameters related to batch 
def set_reference_batch_parameters(config):
  for ref in config['references'].keys():
    try:
      config["references"][ref]['batch']
    except:
      config["references"][ref]['batch'] = None

# set parameters related to benchmarking directory 
def set_benchmark_directory(config,mode):
  import sys 
  for ref in config['references'].keys():
    try:
      config["references"][ref]['output_dir_benchmark']
    except:
      if mode == 'annotation':
        if (config["consensus"]["type"]["CAWPE"]["alpha"][0] != 0) & (config["consensus"]["type"]["CAWPE"]["mode"] != ""):
          sys.exit("@ For running CAWPE the output directory of the benchmarking pipeline should be specified")
        else:
          config["references"][ref]['output_dir_benchmark'] = ""
      else:
        sys.exit("@ In the benchmarking pipeline, all the directories should be specified")
# set parameters related to downsampling 
def set_downsampling_parameters(config):
  for ref in config['references'].keys():
    try: 
      config["references"][ref]['min_cells_per_cluster']
    except: 
      config["references"][ref]['min_cells_per_cluster'] = 0

    try: 
      config["references"][ref]['downsample']['value']
    except: 
      config["references"][ref]['downsample'] = {}
      config["references"][ref]['downsample']['value'] = 0

    try:
      config["references"][ref]['downsample']['stratified']
    except:	
      config["references"][ref]['downsample']['stratified'] = False


# set parameters related to ontology 
def set_ontology_parameters(config,mode):
  import pandas as pd
  
  for ref in config['references'].keys():
    
    try:
      config["references"][ref]["ontology"]["ontology_path"]
      
      try:
        config["references"][ref]["ontology"]["ontology_column"]

      except:
        columns = pd.read_csv(config["references"][ref]["ontology"]["ontology_path"], nrows=0).columns.tolist()
        config["references"][ref]["ontology"]["ontology_column"] = columns[1:]

    except: 
      config["references"][ref]["ontology"] = {}
      if mode == 'annotation':
        config["references"][ref]["ontology"]["ontology_path"] = config["output_dir"] + "/model/" + ref + "/ontology/ontology.csv"
      else:
        config["references"][ref]["ontology"]["ontology_path"] = config["references"][ref]['output_dir_benchmark'] + "/" + ref + '/ontology/ontology.csv'
      config["references"][ref]["ontology"]["ontology_column"] = ['label']
      continue 
    
    if not isinstance(config["references"][ref]["ontology"]["ontology_column"], list): 
      config["references"][ref]["ontology"]["ontology_column"] = [config["references"][ref]["ontology"]["ontology_column"]]
    
    config["references"][ref]["ontology"]["ontology_column"] = ['label'] + config["references"][ref]["ontology"]["ontology_column"]

# return consensus methods 
def get_consensus_methods(config): 
  import sys 
  
  consensus_run = []
  if config["consensus"]["type"]["majority"]["min_agree"][0] != 0:
    consensus_run.append("majority")
  
  if (config["consensus"]["type"]["CAWPE"]["alpha"][0] != 0) & (config["consensus"]["type"]["CAWPE"]["mode"] != ""):
    consensus_run.append("CAWPE")
  
  if len(consensus_run) == 0:
    sys.exit("@ At least one consensus type should be specified")

  return consensus_run

# return the tools to run adding the models to scPred
def get_tools_to_run(config):
  tools_run = config['tools_to_run']
  if "scPred" in tools_run:
    method = config['scPred']['classifier']
    if not isinstance(method, list):
        #Convert into a list
        method = [method]
    tools_to_run = [tool + "_" + m if tool == "scPred" else tool for tool in tools_run for m in (method if tool == "scPred" else [""])]
  else:
    tools_to_run = tools_run
  return(tools_to_run)

# return the tools to run adding the models to scPred
def get_consensus_tools(config):
  consensus_to_run = config['consensus']['tools']
  if 'all' != consensus_to_run:
    if "scPred" in consensus_to_run:
      method = config['scPred']['classifier']
      if not isinstance(method, list):
          #Convert into a list
          method = [method]
      consensus_to_run = [tool + "_" + m if tool == "scPred" else tool for tool in consensus_to_run for m in (method if tool == "scPred" else [""])]
  return(consensus_to_run)

