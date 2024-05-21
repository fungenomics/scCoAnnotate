
# set parameters related to converting gene names 
def set_gene_conversion_parameters(config):
  for ref in config['references'].keys():
    try:
      config["references"][ref]['convert_ref_mm_to_hg']
    except:
      config["references"][ref]['convert_ref_mm_to_hg'] = False


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
      config["references"][ref]['downsample']['per_class']
    except:	
      config["references"][ref]['downsample']['per_class'] = False

# set parameters related to ontology 
def set_ontology_parameters(config):
  import pandas as pd
  
  for ref in config['references'].keys():
    
    try:
      config["references"][ref]["ontology"]["ontology_path"]
      
      try:
        config["references"][ref]["ontology"]["ontology_column"]
      
      except:
        columns = pd.read_csv(config["references"][ref]["ontology"]["ontology_path"], nrows=0).columns.tolist()
        config["references"][ref]["ontology"]["ontology_column"] = columns
        return
    
    except: 
      config["references"][ref]["ontology"] = {}
      config["references"][ref]["ontology"]["ontology_path"] = config['output_dir'] + "/model/" + ref + "/ontology/ontology.csv"
      config["references"][ref]["ontology"]["ontology_column"] = ['label']
      return
    
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


