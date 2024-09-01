#HALLMARK PATHWAY ANALYSIS
{
  HM = gmtPathways('/work/cwt/NACT_OV/1_clustering/genelist/h.all.v2023.2.Hs.symbols.gmt')
  HM
  
  require(GSVA)
  expr = as.matrix(epi_post@assays$RNA@data)
  expr_gsva = gsva(expr,HM,
                   method = 'zscore',parallel.sz=1,kcdf="Gaussian")
  meta = epi_post@meta.data[,c("orig.ident","dcelltype_post")] #类别

  meta <- meta %>% arrange(meta$dcelltype_post)
  data <- expr_gsva[,rownames(meta)]
  
  
  group <- factor(meta[,"dcelltype_post"],ordered = F)
  data1 <-NULL
  dcelltype_post = levels(epi_post$dcelltype_post)
  for(i in dcelltype_post){
    ind <-which(group==i)
    dat <- apply(data[,ind], 1, mean) ##for only stemness, data[ind],others data[,ind]
    data1 <-cbind(data1,dat)
  }
  colnames(data1) <- dcelltype_post
  result<- t(scale(t(data1)))
  
  new_row_names <- gsub("HALLMARK_", "", rownames(result))
  rownames(result) <- new_row_names
  
  
  selected_hallmark = arrange(selected_hallmark,Category)
  selected_hallmark$Name = gsub('REACTIVE_OXIGEN_SPECIES_PATHWAY','REACTIVE_OXYGEN_SPECIES_PATHWAY',  selected_hallmark$Name)
  x = intersect(selected_hallmark$Name,rownames(result))
  result = as.data.frame(result)
  result_selected = result[selected_hallmark$Name,]

  
  category_annotation <- selected_hallmark %>%
    select(Category) 
  
  category_colors <- c("#C8ACD6", "#987D9A", "#0F67B1", "#40534C")
  names(category_colors) <- unique(category_annotation$Category)
  
  # Create the annotation
  left_annotation <- rowAnnotation(
    df = data.frame(Category = category_annotation$Category),
    col = list(Category = category_colors)
  )
  
  dcelltype_post = levels(epi_post@active.ident)
  color.post = c('#C5D8A4','#6FB2D2','#F06161','#3E4A61')
  names(color.post) = dcelltype_post
  top_annotation = HeatmapAnnotation(Cluster = dcelltype_post,
                                     col = list(Cluster = color.post))
  
  
  cell_size = unit(5, "mm")
  ht = Heatmap(result_selected, 
               col = colorRamp2(c(-2, 0, 2), c("#1679AB", "white", "#C80036")),
               width = cell_size * ncol(result_selected),
               height = cell_size * nrow(result_selected),
               show_row_dend = F,
               show_column_dend = F,name = "Z Score",
               cluster_columns = F,cluster_rows = F,
               column_title= 'HALLMARK',
               left_annotation = left_annotation,
               top_annotation = top_annotation,
               show_column_names = T,
               row_names_gp = gpar(fontsize = 8))     
  draw(ht,heatmap_legend_side = "left", annotation_legend_side = "left")
  getwd()
  setwd('/work/cwt/NACT_OV/')
  pdf(file = '/work/cwt/NACT_OV/2_function/epi_post_hallmark.pdf',width = 10,height = 15)
  draw(ht,heatmap_legend_side = "left", annotation_legend_side = "left")
  dev.off()
  
}