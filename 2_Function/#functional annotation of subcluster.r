#functional annotation of subcluster

{
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  
  DE = dplyr::filter(epi.marker2,avg_log2FC > 0.585 & p_val_adj <= 0.05)
  
  TC_DE = lapply(c('Epi0','Epi1','Epi2','Epi3','Epi4'),function(celltype){
    dt = dplyr::filter(DE,cluster == celltype)
    return(dt)
  })
  
  ego1 = enrichGO(gene = rownames(TC_DE[[1]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Hs.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego2 = enrichGO(gene = rownames(TC_DE[[2]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Hs.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego3 = enrichGO(gene = rownames(TC_DE[[3]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Hs.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego4 = enrichGO(gene = rownames(TC_DE[[4]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Hs.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego5 = enrichGO(gene = rownames(TC_DE[[5]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Hs.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )

  plotlist = list()
  
  p1 = dotplot(ego1)+ggtitle('Epi0')
  p2 = dotplot(ego2) +ggtitle('Epi1');p2
  p3 = dotplot(ego3)+ggtitle('Epi2')
  p4 = dotplot(ego4) +ggtitle('Epi3')
  p5 = dotplot(ego5)+ggtitle('Epi4')
  cowplot::plot_grid(p1,p2,p3,p4,p5,ncol = 5)
  
  extractego = function(ego,name){
    sub = ego@result
    sub$Cluster = name
    return(sub)
  }
  
  ego1_df = extractego(ego1,name = 'Epi0')
  ego2_df = extractego(ego2,name = 'Epi1')
  ego3_df = extractego(ego3,name = 'Epi2')
  ego4_df = extractego(ego4,name = 'Epi3')
  ego5_df = extractego(ego5,name = 'Epi4')
  
  
  ego_result = dplyr::bind_rows(ego1_df,ego2_df,ego3_df,ego4_df,ego5_df)
  write.csv(ego_result,'ego_bp_epi_result.csv')
  
  top5_per_cluster <- ego_result %>%
    group_by(Cluster) %>%
    arrange(desc(Count)) %>%
    slice_head(n = 5)
  
  top5_per_cluster$P1 = -log10(top5_per_cluster$p.adjust)
  
  p_epi0 = ggplot(subset(top5_per_cluster,Cluster == 'Epi0'), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#E0AED0") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'Epi0',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  p_epi1 = ggplot(subset(top5_per_cluster,Cluster == 'Epi1'), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#CAEDFF") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'Epi1',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), color = "black", size = 4,hjust = 0.7)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0));p_epi1
  
  p_epi2 = ggplot(subset(top5_per_cluster,Cluster == 'Epi2'), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#FBF0B2") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'Epi2',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), color = "black", size = 4,hjust = 1)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0));p_epi2
  
  p_epi3 = ggplot(subset(top5_per_cluster,Cluster == 'Epi3'), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#FCBAAD") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'Epi3',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), color = "black", size = 4,hjust = 1)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0));p_epi3

  p_epi4 = ggplot(subset(top5_per_cluster,Cluster == 'Epi4'), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#73BBA3") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'Epi4',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), color = "black", size = 4,hjust = 1)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0));p_epi4
  
  cowplot::plot_grid(p_epi0,p_epi1,p_epi2,p_epi3,p_epi4,align = "h",ncol = 3)
  
  ggsave('epi_go_bp.pdf',width= 9.24,height = 5)
  }