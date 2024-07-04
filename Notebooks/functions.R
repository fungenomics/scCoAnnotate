
#--------- HELPER FUNCTIONS ----------------

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
#-------------------------------------------


#------------ THEMES -----------------------

umap_theme = theme(aspect.ratio = 1,
        text = element_text(size = 10), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

default_theme = theme(text = element_text(size = 10), 
                       panel.border = element_rect(colour = "grey90", fill=NA, size=0.5),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.line = element_blank(), 
                       legend.position = "none",
                       strip.background = element_blank())

#-------------------------------------------


#-------- OTHER FUNCTIONS ------------------

get_pred = function(pred, tool, true){
  pred %>%
     select(tool) %>%
     mutate(label = .data[[tool]],
            label = ifelse(!label %in% true$label, NA, label),
            label = factor(label, ordered = TRUE)) %>%
  return()
}

harmonize_unresolved = function(pred, ref_labels){
  pred %>%
  column_to_rownames('cellname') %>%
  mutate(across(where(is.character), ~ifelse(. %in% c(ref_labels$label, 'No Consensus'), ., 'Unresolved'))) %>%
  return()
}

create_color_pal = function(class, mb = 'Juarez'){
  pal = sample(met.brewer(mb, length(class)))  
  names(pal) = class
  pal['Unresolved'] = 'lightgrey'
  pal['No Consensus'] = 'grey'
  return(pal)
}

calculate_percentage_unresolved = function(pred, order){
  warn = pred %>% 
    select(order) %>%
    pivot_longer(order) %>% 
    mutate(value = factor(value)) %>%
    group_by(name) %>%
    count(value, .drop = F) %>%
    mutate(frac = n/sum(n)*100) %>%
    filter(!(!name == 'Consensus'  & value == 'No Consensus')) %>%
    filter(value %in% c('No Consensus', 'Unresolved')) %>%
    mutate(warn = case_when(frac >= 70 ~ 'HIGH', 
                            frac < 70 & frac > 30 ~ 'MEDIUM',
                            frac <= 30 ~ 'LOW'))
   if(nrow(warn) != 0){
   warn = data.frame(TOOL = warn$name, 
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


#------------ PLOTS -------------------------

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

# Plot class stat per fold (F1 etc) as barplot
plot_stat = function(cm_byclass, stat){

p = cm_byclass %>% 
  as.data.frame() %>%
  rownames_to_column('class') %>%
  separate(class, into = c(NA, 'class'), sep = ': ') %>%
  ggplot(aes(reorder(class, -.data[[stat]]), .data[[stat]])) +
  geom_bar(stat = 'identity', col = 'white', fill = 'lightgrey') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        aspect.ratio = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = c(1, 0.5), linetype = 'dotted', color = 'red') +
  default_theme

return(p)
}

# plot F1 accross folds for each class as a boxplot  
plot_stat_boxplot = function(list, tool, stat){
  
df = lapply(list[[tool]], get_stat, stat = stat) %>% bind_rows()

df[is.na(df)] = 0

df %>%
  ggplot(aes(reorder(class, -.data[[stat]], mean), .data[[stat]])) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        aspect.ratio = 0.5) +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0)) +
  geom_hline(yintercept = c(1, 0.5), linetype = 'dotted', color = 'red') +
  default_theme
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

split = c('Consensus', rep('tools', length(tools)-1))

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
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1)) +
  default_theme

  return(b)
}

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
  theme(axis.title = element_blank(),
        aspect.ratio = 0.5) +
  default_theme

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

#-------------------------------------------