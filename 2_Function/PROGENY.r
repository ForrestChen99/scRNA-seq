##PROGENY
#240831 @Forrest Chen
#Dsecribes the Progeny workflow of single cell sequenccing
{
  #devtools::install_github("saezlab/progeny")
  library(progeny)
  
  
  CellsClusters <- data.frame(Cell = names(Idents(tumor_ident)), 
                              CellType = as.character(Idents(tumor_ident)),
                              stringsAsFactors = FALSE)
  
  progeny_obj = progeny(tumor_ident,scale = F,organism = 'Mouse',top = 500,return_assay = TRUE)
  
  progeny_obj <- Seurat::ScaleData(progeny_obj, assay = "progeny") 
  
  progeny_scores_df <- 
    as.data.frame(t(GetAssayData(progeny_obj, slot = "scale.data", 
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 
  
  
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters,by = 'Cell')
  
  summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
  
  paletteLength = 100
  myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
  progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(summarized_progeny_scores_df)/paletteLength, 
                        max(summarized_progeny_scores_df), 
                        length.out=floor(paletteLength/2)))
  
  pheatmap(summarized_progeny_scores_df[,-1],fontsize=14, 
           fontsize_row = 10, 
           main = "PROGENy",color =colorRampPalette(c("Darkblue", "white","#E64B35FF"))(100),
           breaks = progenyBreaks,cluster_rows = F,
           treeheight_col = 0,  border_color = NA)
  
  write_rds(summarized_progeny_scores_df,'Progeny.rds.gz')
  #summarized_progeny_scores_df = read_rds('Progeny.rds.gz')
  
  result1 = as.data.frame(t(result))
  result2 = cbind(summarized_progeny_scores_df,result1)
  paletteLength = 100
  breaks2 = c(seq(min(result2), 0, 
                  length.out=ceiling(paletteLength/2) + 1),
              seq(max(result2)/paletteLength, 
                  max(result2), 
                  length.out=floor(paletteLength/2)))
  result3 = result2[,-c((ncol(result2) - 2):ncol(result2))]
  
  pdf('PROGENy_defined.heatmap.pdf',width = 8,height = 3)
  ComplexHeatmap::pheatmap(result3,
                           cluster_rows = F,
                           cluster_cols = F,
                           show_rownames = T,
                           show_colnames = T,
                           scale = 'column',
                           color =colorRampPalette(c("Darkblue", "white","#E64B35FF"))(100),breaks = breaks2,
                           #cellwidth = 15, cellheight = 15,
                           #fontsize = 10
  )
  dev.off()
  
}