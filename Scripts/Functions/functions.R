

#-------- THEME --------------------------------------------

umap_theme = theme(
        aspect.ratio = 1,
        text = element_text(size = 10), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        panel.border = element_rect(colour = "grey90", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
        )

notebook_theme = theme(
        text = element_text(size = 10), 
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        panel.border = element_rect(colour = "grey90", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        )

#-------- PREPROCESSING ------------------------------------
# load and rename the Rdata objects
loadRDa <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# get the expression matrix and labels from the reference. 
# check if it's a .csv, Seurat or SingleCellExperiment object.
get_data_reference <- function(ref_path,
                               lab_path,
                               batch_path = NULL){
  ref_ext <- tools::file_ext(ref_path) %>% str_to_lower
  if(ref_ext == 'csv'){ #If the expression is csv it assumes that the labels is csv
    if(tools::file_ext(lab_path) != 'csv'){
      stop('@ Labels should be provided as a .csv data.frame')
    }
    
    # read reference 
    mtx = data.table::fread(ref_path, header = T) %>% column_to_rownames('V1')
    #read labels   
    lab = data.table::fread(lab_path, header = T) %>% column_to_rownames('V1')  
  } else if(ref_ext %in% c('rda','rdata','rds')){ #If the rda rdata rds it assumes that the label is a vector of column name
    if(ref_ext %in% c('rda','rdata')){
      sce <- loadRDa(ref_path)
    } else if(ref_ext %in% c('rds')){
      sce <- readRDS(ref_path)
    }
    if(is(sce,'Seurat')){
      mtx <- sce@assays$RNA@counts %>% as.matrix %>% t() %>% as.data.frame()
      lab <- data.frame(row.names = colnames(sce),
                        label = sce@meta.data[,lab_path,drop=T])
      
      if(!is.null(batch_path)){
        lab$batch <- sce@meta.data[,batch_path,drop=T]
      }
      
    } else if(class(sce) %in% c('SingleCellExperiment','LoomCellExperiment')){
      mtx <- assay(sce,'counts') %>% as.matrix() %>% t() %>% as.data.frame()
      lab <- data.frame(row.names = colnames(sce),
                        label = colData(sce)[,lab_path,drop=T])
      
      if(!is.null(batch_path)){
        lab$batch <- colData(sce)[,batch_path,drop=T]
      }
      
    } else{
      stop("@ Object is not Seurat (v4) nor SingleCellExperiment")
    }
  } else{
    stop('@ Formats provided are not compatible. Acceptable formats are .csv, Seurat (v4) object, SingleCellExperiment object (.rds,.rdata)')
  }
  return(list(exp = mtx,
              lab = lab)
  )
}

# get the expression matrix from the query and extract them if the input is a Seurat or SingleCellExperiment object
get_data_query <- function(query_path){
  query_ext <- tools::file_ext(query_path) %>% str_to_lower
  if(query_ext == 'csv'){ #If the expression is csv it assumes that the labels is csv
    # read reference 
    mtx = data.table::fread(query_path, header = T) %>% column_to_rownames('V1')
  } else if(query_ext %in% c('rda','rdata','rds')){ #If the rda rdata rds it assumes that the label is a vector of column name
    if(query_ext %in% c('rda','rdata')){
      sce <- loadRDa(query_path)
    } else if(query_ext %in% c('rds')){
      sce <- readRDS(query_path)
    }
    if(is(sce,'Seurat')){
      mtx <- sce@assays$RNA@counts %>% as.matrix %>% t %>% as.data.frame()
    } else if(class(sce) %in% c('SingleCellExperiment','LoomCellExperiment')){
      mtx <- assay(sce,'counts') %>% as.matrix() %>% t() %>% as.data.frame()
    } else{
      stop("@ Object is not Seurat (v4) nor SingleCellExperiment")
    }
  } else{
    stop('@ Formats provided are not compatible. Acceptable formats are .csv, Seurat object, SingleCellExperiment object (.rds,.rdata)')
  }
  return(mtx)
}

# convert mouse to human gene names 
mapfun = function(mousegenes){
    gns    = mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
    mapped = AnnotationDbi::select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
    naind  = is.na(mapped$Homo_sapiens)
    hsymb  = mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
    out    = data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
    out$Human_symbol[!naind] = hsymb
    return(out)
}

downsample = function(labels, downsample_stratified, downsample_value = 1){
  
  downsample_stratified = if(downsample_stratified) "label" else NULL
  
  if(downsample_value >= 1){
      labels = labels %>%  
               rownames_to_column('cell') %>%
               group_by(across(all_of(downsample_stratified))) %>% 
               dplyr::slice_sample(n = downsample_value,replace = F) %>%
               column_to_rownames('cell')
  }else{
      labels = labels %>% 
               rownames_to_column('cell') %>%
               group_by(across(all_of(downsample_stratified))) %>% 
               dplyr::slice_sample(prop = downsample_value,replace = F) %>%
               column_to_rownames('cell')
  }

  return(labels)
}

remove_small_clusters = function(labels, min_cells){

  rmv_labels = names(which(table(labels$label) < min_cells))
  labels = labels %>% filter((!label %in% rmv_labels))
  
  message(paste0(paste0(rmv_labels,collapse = '-'),
  	' classes were remove because of lower number of cells (< ',as.character(min_cells),')'))

  return(labels)
}

#-------- CONSENSUS ------------------------------------

read_prediction_files = function(path, tools){
	
	files = list.files(path, pattern = 'pred.csv', recursive = T, full.names = T)

    l = list()
    
    for(f in files){
    	l[[basename(dirname(f))]] = data.table::fread(f)
    }

    consensus = l %>% 
          reduce(left_join, by = "cell") %>% 
          rename('cellname' = 'cell') %>% 
          select(all_of(c('cellname', tools))) %>% 
          as.data.frame()

    return(consensus)
}


## This function takes the data.frame ontology and a vector (pred)
## and convert the original label (from) to its ontology (to). 
## The output is the vector
apply_ontology <- function(df_ontology,
                           pred,
                           from = "label",
                           to){
  df_ontology = as.data.frame(df_ontology)
  ont <- setNames(nm = as.character(df_ontology[,from,drop=T]), object = as.character(df_ontology[,to,drop=T])
  )

  match_ont <- as.character(ont[pred])
  pred[!is.na(match_ont)] <- match_ont[!is.na(match_ont)]
  
  return(pred)
}

harmonize_unresolved = function(pred, ref_labels){
  pred %>%
  column_to_rownames('cellname') %>%
  mutate(across(where(is.character), ~ifelse(. %in% c(ref_labels), ., 'Unresolved'))) %>%
  rownames_to_column('cellname') %>%
  return()
}

get_consensus = function(v, min = 2){
  
  uniqv = unique(v)
  matches = tabulate(match(v, uniqv))
  max_match = max(matches)
  
  ties = ifelse(length(which(matches == max_match)) > 1, T, F)
  
  if (max_match < min) {
    return("No Consensus")
  } else if (ties) {
    return("No Consensus") 
  } else {
    return(uniqv[which.max(matches)])
  }
}

get_max = function(v){
  uniqv = unique(v)
  matches = tabulate(match(v, uniqv))
  max_match = max(matches) / sum(matches)
  return(max_match)
}

get_entropy = function(v){
  uniqv = unique(v)
  matches = tabulate(match(v, uniqv))
  matches <- matches / sum(matches)
  entropy = -sum(matches*log2(matches),na.rm = T)
  ## Normalizing by the number of labels to make it  [0,1]
  entropy = entropy / log2(length(matches))
  #If NAN is because it assgined everything to the same category
  if(is.na(entropy)){
    entropy <- 0
  }
  return(entropy)
}

CAWPE = function(x, alpha = 4){
  (as.numeric(x["mean_metric"])^alpha)*as.numeric(x['prob'])
}


get_cawpe_columns = function(CAWPE_type){
    if(CW_tp == 'CAWPE_T'){
      cols = c('tool')
    }else if(CW_tp == 'CAWPE_CT'){
      cols = c('tool', 'class')
    }
  return(cols) 
}

check_ontology_hierarchy = function(ontology){ 
  if(!length(unique(ontology$label)) == nrow(ontology)){
    stop(paste0("Ontology does not have correct hierarchy"))
  }
}

#----- NOTEBOOK FUNCTIONS -----------------------------------

create_color_pal = function(class, mb = 'Juarez'){
  pal = sample(met.brewer(mb, length(class)))  
  names(pal) = class
  pal['Unresolved'] = 'lightgrey'
  pal['No Consensus'] = 'grey'
  return(pal)
}

calculate_percentage_unresolved = function(pred, order, cons_tools){
  warn = pred %>% 
    select(order) %>%
    pivot_longer(order) %>% 
    mutate(value = factor(value)) %>%
    group_by(name) %>%
    count(value, .drop = F) %>%
    mutate(frac = n/sum(n)*100) %>%
    filter(!(!name == 'Consensus' & value == 'No Consensus')) %>%
    filter(value %in% c('No Consensus', 'Unresolved')) %>%
    mutate(in_cons = ifelse(name %in% cons_tools & name != 'Consensus', 'YES', '')) %>%
    mutate(warn = case_when(frac >= 70 ~ 'HIGH', 
                            frac < 70 & frac > 30 ~ 'MEDIUM',
                            frac <= 30 ~ 'LOW'))
   if(nrow(warn) != 0){
   warn = data.frame(TOOL = warn$name, 
                     'IN CONSENSUS' = warn$in_cons,
                     LABEL = warn$value, 
                     PERCENTAGE = warn$frac, 
                     FLAG = warn$warn) %>%
          mutate(TOOL = factor(TOOL, levels = order)) 

   warn = warn[order(warn$TOOL),]

   warn$FLAG = cell_spec(warn$FLAG, 
                         bold = T,
                         background = case_when(warn$FLAG == 'HIGH' ~ "red",
                                                warn$FLAG == 'MEDIUM' ~ "yellow",
                                                warn$FLAG == 'LOW' ~ "green"))
   
   warn$TOOL = cell_spec(warn$TOOL,
                         bold = ifelse(warn$TOOL == 'Consensus', T, F),
                         background = ifelse(warn$TOOL == 'Consensus', 'black', 'white'),
                         color = ifelse(warn$TOOL == 'Consensus', 'white', 'black'))
   }
   
   return(warn)
}


get_pred = function(pred, tool, true){
  pred %>%
     select(tool) %>%
     mutate(label = .data[[tool]],
            label = ifelse(!label %in% true$label, NA, label),
            label = factor(label, ordered = TRUE)) %>%
  return()
}

# gets stat for each fold and returns data frame 
get_stat = function(x, stat){
  x$byClass %>% 
  as.data.frame() %>%
  rownames_to_column('class') %>%
  separate(class, into = c(NA, 'class'), sep = ': ') %>%
  select(class, .data[[stat]]) %>%
  mutate(fold = x$fold,
         tool = x$tool)
}

# gets all stats for each fold and returns data frame 
get_all_stats = function(x){
  x$byClass %>% 
  as.data.frame() %>%
  rownames_to_column('class') %>%
  separate(class, into = c(NA, 'class'), sep = ': ') %>%
  mutate(fold = x$fold,
         tool = x$tool)
}

#----- PLOTS FOR NOTEBOOK ANNOTATION ------------------------------------

plot_tool_correlation_heatmap = function(seurat, tools){
 
 mat = query@meta.data %>%
  select(all_of(tools)) %>%
  rownames_to_column('cell') %>%
  pivot_longer(!cell) %>%
  mutate(value = factor(value)) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames('cell') %>%
  cor()

  mat[is.na(mat)] = 0

  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("#2274A5", "beige", "#F75C03"))

 count = query@meta.data %>%
  select(all_of(tools)) %>%
  rownames_to_column('cell') %>%
  pivot_longer(!cell) %>%
  filter(value %in% c('Unresolved', 'No Consensus')) %>%
  dplyr::count(name, .drop = F) %>%
  mutate(freq = round(n/nrow(query@meta.data)*100)) %>%
  select(!n) %>%
  column_to_rownames('name')

  count[setdiff(names(seurat@meta.data %>% select(tools)), rownames(count)),] = 0
  count = count[order(match(rownames(count), colnames(mat))), , drop = FALSE]

  ha = columnAnnotation('% Unresolved/No Consensus' = anno_barplot(count, border = F, gp = gpar(fill = '#596475', col = '#596475')))

  h = ComplexHeatmap::Heatmap(mat,
                        name = 'Correlation',
                        col = col_fun,
                        width = ncol(mat)*unit(7, "mm"),
                        height = nrow(mat)*unit(7, "mm"),
                        rect_gp = gpar(col = "white", lwd = 2), 
                        top_annotation = ha, 
                        show_column_dend = F)
  return(h)
}

plot_bar_largest_group = function(seurat, meta_column = '', pal = pal, fr = 0.1){
  
df = seurat@meta.data %>%
  count(seurat_clusters, .data[[meta_column]]) %>%
  group_by(seurat_clusters) %>% 
  mutate(`%` = (n / sum(n)))  %>% 
  mutate(meta = ifelse(`%` < fr, NA, .data[[meta_column]]))

pal = pal[unique(na.omit(df$meta))]

df = df %>% 
     mutate(meta = factor(meta, levels = c(NA, names(pal)), exclude = NULL))

p1 = df %>%
  ggplot(aes(x = seurat_clusters, y = `%`, fill = meta, text = sprintf(" %s <br> %s ", 
                                                                       meta, 
                                                                       scales::percent(`%`, scale = 100, accuracy = .1)))) +
  geom_bar(stat = 'identity', position="fill") +
  scale_fill_manual(values = pal, na.value = 'white', name = '') +
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format(scale = 100, accuracy = 1)) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.title = element_blank(),
        axis.line = element_line(size = 0.5),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        aspect.ratio = 0.5) +
  notebook_theme

 p2 = ggplotly(p1, tooltip = c('text')) %>% toWebGL()
 
 return(p2)
}

plot_percentage_predicted_consensus_class = function(seurat, tools){
  col_fun = circlize::colorRamp2(c(0, 50, 100), c("#2274A5", "beige", "#F75C03"))

x = seurat@meta.data %>%
  select(c('Consensus', tools)) %>%
  group_by(Consensus) %>% 
  pivot_longer(!c('Consensus')) %>%
  group_by(Consensus, name) %>%
  count(value, .drop = F) %>% 
  mutate(freq = n/sum(n)*100) %>%
  filter(Consensus == value) %>%
  select(name, freq, Consensus) %>%
  pivot_wider(values_from = 'freq', names_from = 'Consensus',values_fill = 0) %>%
  column_to_rownames('name')  

  h = ComplexHeatmap::Heatmap(x,
                        name = '%',
                        col = col_fun,
                        width = ncol(x)*unit(4, "mm"),
                        height = nrow(x)*unit(4, "mm"),
                        rect_gp = gpar(col = "white", lwd = 2), row_names_side = 'left',
                        show_row_dend = F, column_names_gp = gpar(size = 7))
  
  return(h)
}

color_class_seurat = function(seurat, meta_column, pal){
  list = list()
  pal['Unresolved'] = 'red'
  pal['No Consensus'] = 'red'
  Idents(seurat) = meta_column
  class = (table(query@meta.data[[meta_column]]) %>% as.data.frame() %>% filter(Freq > 20))$Var
  
  for(c in class){
    lab =  names(Idents(seurat)[Idents(seurat) == c])
    p = DimPlot(seurat, cells.highlight = lab,  cols = 'lightgrey', cols.highlight = pal[c], pt.size = 1) + umap_theme + ggtitle(c)
    list[[c]] = p 
  }
  
  return(list)
}

feature_plot_seurat = function(seurat, genes){
  list = list()

  genes = genes[genes %in% rownames(seurat@assays$RNA)]
  
  for(g in genes){
    p = FeaturePlot(seurat, features = g, cols = c("#F2EFC7", "#BC412B"), order = T) + umap_theme + ggtitle(g)
    list[[g]] = p 
  }  
  return(list)
}

feature_plot_seurat_meta = function(seurat, meta_cols){
  list = list()

  meta_cols = meta_cols[meta_cols %in% colnames(seurat@meta.data)]
  for(mc in meta_cols){
    p = FeaturePlot(seurat, features = mc, cols = c("#F2EFC7", "#BC412B"), order = T) + umap_theme + theme(legend.position = 'right') + ggtitle(mc)
    list[[mc]] = p 
  }  
  return(list)
}

umap_plotly = function(seurat, meta_column, pal){

  p1 = cbind(seurat@reductions$umap@cell.embeddings, seurat@meta.data) %>%
  slice(sample(1:n())) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = .data[[meta_column]], text = .data[[meta_column]])) + 
  geom_point(alpha = 0.8) + 
  scale_color_manual(values = pal) +
  theme_bw() + umap_theme + theme(legend.position = 'right') 

  p2 = ggplotly(plot = p1, tooltip = c('text')) %>% layout(autosize = F, width = 550, height = 450) %>% toWebGL()

  return(p2)
}


#----- PLOTS FOR NOTEBOOK BENCHMARK ------------------------------------

# Plot confusion matrix as a heatmap 
plot_cm = function(cm_table){
  col_fun = circlize::colorRamp2(c(range(cm_table)[1], 
                                   range(cm_table)[2]/2, 
                                   range(cm_table)[2]), 
                                 c("#5C80BC", "#F2EFC7", "#FF595E")) 
  
  h = Heatmap(cm_table,
              name = 'Counts',
              col = col_fun,
              width = ncol(cm_table)*unit(2, "mm"),
              height = nrow(cm_table)*unit(2, "mm"),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 7), 
              column_title = 'True Class', 
              row_title = 'Predicted Class')
  
  return(h)
}

# Plot class stat per fold (metric) as barplot
plot_stat = function(cm_byclass, stat){

p = cm_byclass %>% 
  as.data.frame() %>%
  rownames_to_column('class') %>%
  separate(class, into = c(NA, 'class'), sep = ': ') %>%
  ggplot(aes(reorder(class, -.data[[stat]]), .data[[stat]])) +
  geom_bar(stat = 'identity', col = 'white', fill = 'lightgrey') +
  theme_bw() +
  theme(text = element_text(size = 10), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.5), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        aspect.ratio = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = c(1, 0.5), linetype = 'dotted', color = 'red')  +
  notebook_theme

return(p)
}

# plot metric accross folds for each class as a boxplot  
plot_stat_boxplot = function(list, tool, stat){
  
df = lapply(list[[tool]], get_stat, stat = stat) %>% bind_rows()

df[is.na(df)] = 0

df %>%
  ggplot(aes(reorder(class, -.data[[stat]], mean), .data[[stat]])) +
  geom_boxplot() +
  theme_bw() +
  theme(text = element_text(size = 10), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.5), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        aspect.ratio = 0.5) +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0)) +
  geom_hline(yintercept = c(1, 0.5), linetype = 'dotted', color = 'red')  +
  notebook_theme
}

# plot average stat for all tools 
plot_mean_tool = function(list, stat, tools){

df = lapply(list, function(x){lapply(x, get_stat, stat = stat) %>% bind_rows()})

df = bind_rows(df) %>% 
  group_by(class, tool) %>%
  mutate(mean = mean(.data[[stat]])) %>%
  distinct(class, tool, mean) %>% 
  pivot_wider(names_from = 'class', values_from = mean) %>%
  column_to_rownames('tool')

df[is.na(df)] = 0

col_fun = circlize::colorRamp2(c(0, 
                                range(df)[2]/2, 
                                range(df)[2]), 
                                 c("#3B5B91", "#F2EFC7", "#CC0007")) 
cons_number = length(grep(pattern = "^Consensus_",x = tools))
split = c(rep("Consensus",cons_number), rep('tools', length(tools)-cons_number))

h = Heatmap(df,
            name = paste('Mean ', stat),
            col = col_fun,
            width = ncol(df)*unit(4, "mm"),
            height = nrow(df)*unit(6, "mm"),
            row_names_side = 'left',
            row_names_gp = gpar(fontsize = 12),
            show_column_dend = F,
            show_row_dend = F, 
            row_split = split,
            cluster_row_slices = F, 
            row_title = NULL)

return(h)
}

plot_n_cells_per_class = function(df){
 mean_n = df %>%
  count(label, fold) %>%
  group_by(label) %>%
  summarise(mean = round(mean(n)))

 b = df %>%
  count(label, fold) %>%
  ggplot(aes(reorder(label, desc(mean)), mean)) +
  geom_bar(data = mean_n, 
           mapping = aes(reorder(label, desc(mean)), mean), 
           stat = 'identity', 
           fill = 'grey90') +
  geom_text(data = mean_n, mapping = aes(label = mean), hjust = -0.2, vjust = -0.2, angle = 45) +
  geom_point(aes(label, n), alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(mean_n$mean) + (max(mean_n$mean)*0.2))) +
  ylab('N') +
  theme(text = element_text(size = 10), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.5), 
        panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1))  +
  notebook_theme

  return(b)
}

#-----------------------------------------------------------------------
