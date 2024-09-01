#cell chat
##Cell chat analysis of CAF and tumor cells
{
  DotPlot(fib,features=c('ANTXR1'))
  epi_fib = merge(epi,fib)
  table(epi_fib$dcelltype)
  
  
  library(cellchat)
  table(epi_fib@active.ident)
  epi_fib_ds = subset(epi_fib,downsample = 300)
  
  data.input <- GetAssayData(epi_fib_ds, assay = "RNA", slot = "data")
  identity <- subset(epi_fib_ds@meta.data, select = "dcelltype")
  cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "dcelltype")
  
  cellchatDB <- cellchatDB.human
  {
    showDatabaseCategory(cellchatDB)
    colnames(cellchatDB$interaction)
    cellchatDB$interaction[1:4,1:4]
    head(cellchatDB$cofactor)
    head(cellchatDB$complex)
    head(cellchatDB$geneInfo)
  }
  cellchatDB.use <- subsetDB(cellchatDB, search = "Secreted Signaling")
  cellchat@DB <- cellchatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  unique(cellchat@idents)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  ##每个细胞如何跟别的细胞互作（互作的强度或概率图）
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_count.pdf',width = 22,height = 13)
  mat <- cellchat@net$count
  par(mfrow = c(3,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
 
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_weight.pdf',width = 22,height = 13)
  mat <- cellchat@net$weight
  par(mfrow = c(3,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  
  getwd()
  write_rds(cellchat,'/work/cwt/NACT_OV/3_cellchat/cellchat_epi_fib.rds.gz',compress = 'gz')
  
  #F8_CLEC3B AND EPI1? SPECIFIC PATHWAY
  df.net <- subsetCommunication(cellchat)
  df.net1 = subsetCommunication(cellchat,sources.use = c(13),targets.use = c(2))
  
  
  p1 = netVisual_bubble(cellchat,sources.use = c(13),targets.use = c(1:5),remove.isolate = FALSE)
  p2 = netVisual_bubble(cellchat,sources.use = c(1:5),targets.use = c(13),remove.isolate = F)
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_bubble_Secret.pdf',width = 9.5,height = 7.4)
  p1+p2+plot_layout(guides = 'collect')
  dev.off()
  
  p3 = netVisual_bubble(cellchat,sources.use = c(6:13),targets.use = c(1:5),remove.isolate = FALSE)
  p4 = netVisual_bubble(cellchat,sources.use = c(1:5),targets.use = c(6:13),remove.isolate = FALSE)
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_bubble_secret_whole.pdf',width = 10.22,height = 10.8)
  p3+p4+plot_layout(guides = 'collect',ncol = 1)
  dev.off()
  
  
  
  
  netVisual_bubble(cellchat,sources.use = c(6:13),targets.use = c(1:5),remove.isolate = FALSE)
  

  
  ```````````````````Using ECM-Receptor for cellchat```````````````````````````````````
  cellchat_ECM = createCellChat(object = data.input, meta = identity,  group.by = "dcelltype")
  
  cellchatDB.ECM <- subsetDB(CellChatDB, search = "ECM-Receptor")
  cellchat_ECM@DB <- cellchatDB.ECM # set the used database in the object
  
  cellchat_ECM <- subsetData(cellchat_ECM)
  cellchat_ECM <- identifyOverExpressedGenes(cellchat_ECM)
  cellchat_ECM <- identifyOverExpressedInteractions(cellchat_ECM)
  cellchat_ECM <- projectData(cellchat_ECM, PPI.human)
  
  unique(cellchat_ECM@idents)
  cellchat_ECM <- computeCommunProb(cellchat_ECM, raw.use = TRUE)
  cellchat_ECM <- filterCommunication(cellchat_ECM, min.cells = 3)
  
  cellchat_ECM <- computeCommunProbPathway(cellchat_ECM)
  cellchat_ECM <- aggregateNet(cellchat_ECM)
  
  groupSize <- as.numeric(table(cellchat_ECM@idents))
  
  netVisual_circle(cellchat_ECM@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_ECM@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  netVisual_bubble(cellchat_ECM,sources.use = c(6:13),targets.use = c(1:5),remove.isolate = FALSE)
  netVisual_bubble(cellchat_ECM,sources.use = c(1:5),targets.use = c(6:13),remove.isolate = FALSE)
  
  write_rds(cellchat_ECM,'/work/cwt/NACT_OV/3_cellchat/cellchat_epi_fib_ECM.rds.gz',compress = 'gz')

  
  p1 = netVisual_bubble(cellchat_ECM,sources.use = c(13),targets.use = c(1:5),remove.isolate = FALSE)
  p2 = netVisual_bubble(cellchat_ECM,sources.use = c(1:5),targets.use = c(13),remove.isolate = F)
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_bubble_ecm.pdf',width = 10.78,height = 11.3)
  p1+p2+plot_layout(guides = 'collect',width = c(1,1))
  dev.off()
  
  p3 = netVisual_bubble(cellchat_ECM,sources.use = c(6:13),targets.use = c(1:5),remove.isolate = FALSE)
  p4 = netVisual_bubble(cellchat_ECM,sources.use = c(1:5),targets.use = c(6:13),remove.isolate = FALSE)
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_bubble_ecm_whole.pdf',width = 10.22,height = 17)
  p3+p4+plot_layout(guides = 'collect',ncol = 1,heights = c(1.5,1))
  dev.off()
  
  
  DotPlot(epi_fib,features = c('NEGR1','ADGRE5','CD55','CD46','JAG1','MPZL1','APP','CD74'))
  
  
  ``````````````````Using CELL-CELL-CONTACT for cellchat```````````````````````````````````
  cellchat_CC = createCellChat(object = data.input, meta = identity,  group.by = "dcelltype")
  
  cellchatDB.CC <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
  cellchat_CC@DB <- cellchatDB.CC # set the used database in the object
  
  cellchat_CC <- subsetData(cellchat_CC)
  cellchat_CC <- identifyOverExpressedGenes(cellchat_CC)
  cellchat_CC <- identifyOverExpressedInteractions(cellchat_CC)
  cellchat_CC <- projectData(cellchat_CC, PPI.human)
  
  unique(cellchat_CC@idents)
  cellchat_CC <- computeCommunProb(cellchat_CC, raw.use = TRUE)
  cellchat_CC <- filterCommunication(cellchat_CC, min.cells = 3)
  
  cellchat_CC <- computeCommunProbPathway(cellchat_CC)
  cellchat_CC <- aggregateNet(cellchat_CC)
  
  groupSize <- as.numeric(table(cellchat_CC@idents))
  
  p3 = netVisual_bubble(cellchat_CC,sources.use = c(6:13),targets.use = c(1:5),remove.isolate = FALSE)
  p4 = netVisual_bubble(cellchat_CC,sources.use = c(1:5),targets.use = c(6:13),remove.isolate = FALSE)
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_bubble_CC_whole.pdf',width = 10,height = 7.7)
  p3+p4+plot_layout(guides = 'collect',ncol = 1)
  dev.off()
  
  
  netVisual_bubble(cellchat_CC,sources.use = c(13,1:5),targets.use = c(1:5,13),remove.isolate = T)
  
  p1 = netVisual_bubble(cellchat_CC,sources.use = c(13),targets.use = c(1:5),remove.isolate = FALSE)
  p2 = netVisual_bubble(cellchat_CC,sources.use = c(1:5),targets.use = c(13),remove.isolate = F)
  
  
  
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_bubble_CC.pdf',width = 9.5,height = 7.4)
  p1+p2+plot_layout(guides = 'collect',width = c(1,1))
  dev.off()
  
}

#CELLCHAT FOR POST AND PRE CHEMOTHERAPY ANALYSIS
{
  library(cellchat)
  #epi_fib_ds = subset(epi_fib,downsample = 300)
  
  epifib_pre = subset(epi_fib_ds,subset = treatment_phase == 'naive')
  pre.data = GetAssayData(epifib_pre, assay = "RNA", slot = "data")
  identity <- subset(epifib_pre@meta.data, select = "dcelltype")
  cellchat_pre <- createCellChat(object = pre.data, meta = identity,  group.by = "dcelltype")
  
  
  epifib_post = subset(epi_fib_ds,subset = treatment_phase == 'post')
  data.input <- GetAssayData(epifib_post, assay = "RNA", slot = "data")
  identity <- subset(epifib_post@meta.data, select = "dcelltype")
  cellchat_post <- createCellChat(object = data.input, meta = identity,  group.by = "dcelltype")
  
  
  cellchatDB <- CellChatDB.human

  cellchat_pre@DB <- cellchatDB # set the used database in the object
  cellchat_post@DB <- cellchatDB
  
  #pre
  cellchat_pre <- subsetData(cellchat_pre)
  cellchat_pre <- identifyOverExpressedGenes(cellchat_pre)
  cellchat_pre <- identifyOverExpressedInteractions(cellchat_pre)
  cellchat_pre <- projectData(cellchat_pre, PPI.human)
  
  unique(cellchat_pre@idents)
  cellchat_pre <- computeCommunProb(cellchat_pre, raw.use = TRUE)
  cellchat_pre <- filterCommunication(cellchat_pre, min.cells = 3)
  
  cellchat_pre <- computeCommunProbPathway(cellchat_pre)
  cellchat_pre <- aggregateNet(cellchat_pre)
  
  
  #post
  cellchat_post <- subsetData(cellchat_post)
  cellchat_post <- identifyOverExpressedGenes(cellchat_post)
  cellchat_post <- identifyOverExpressedInteractions(cellchat_post)
  cellchat_post <- projectData(cellchat_post, PPI.human)
  
  unique(cellchat_post@idents)
  cellchat_post <- computeCommunProb(cellchat_post, raw.use = TRUE)
  cellchat_post <- filterCommunication(cellchat_post, min.cells = 3)
  
  cellchat_post <- computeCommunProbPathway(cellchat_post)
  cellchat_post <- aggregateNet(cellchat_post)
  
  
  #combined
  object.list <- list(Pre = cellchat_pre, Post = cellchat_post)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  #save(object.list, file = "/work/cwt/NACT_OV/3_cellchat/cellchat_object.list_epifib.RData")
  #save(cellchat, file = "/work/cwt/NACT_OV/3_cellchat/cellchat_merged_epifib.RData")
  
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  #> 
  pdf('/work/cwt/NACT_OV/3_cellchat/epi_fib_hm_postvspre.pdf',width = 10.13,height = 4.4)
  gg1 + gg2
  dev.off()
 
  netVisual_bubble(cellchat, sources.use = 13, targets.use = c(1:5),  comparison = c(1, 2), angle.x = 45)
  
  gg1 <- netVisual_bubble(cellchat, sources.use = c(13), targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  gg2 <- netVisual_bubble(cellchat, sources.use = 13, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  gg1 + gg2
  
  
  
  ##Identify dysfunctional signaling by using differential expression anlaysis
  #240823：这里的问题在于她卡得阈值不一定是针对所感兴趣的亚群。所以针对某一亚群，要从net里面挑出来具体的。
  
  {
    # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
    pos.dataset = "Post"
    # define a char name used for storing the results of differential expression analysis
    features.name = paste0(pos.dataset, ".merged")
    
    # perform differential expression analysis 
    # Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 
    
    cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05) 
    #> Use the joint cell labels from the merged CellChat object
    
    # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
    net <- netMappingDEG(cellchat, features.name = features.name)
    # extract the ligand-receptor pairs with upregulated ligands in post
    net.up <- subsetCommunication(cellchat, net = net, datasets = "Post",ligand.logFC = 0.05, receptor.logFC = NULL)
    # extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in pre, i.e.,downregulated in LS
    net.down <- subsetCommunication(cellchat, net = net, datasets = "Pre",ligand.logFC = -0.05, receptor.logFC = NULL)
 
    gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
    gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
    
    #bubble plot visualization
    {
      pairLR.use.up = net.up[, "interaction_name", drop = F]
      gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(6:13), targets.use = c(1:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
      #> Comparing communications on a merged object
      pairLR.use.down = net.down[, "interaction_name", drop = F]
      gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 13, targets.use = c(1:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
      #> Comparing communications on a merged object
      gg1 + gg2
    }
    
    #F5_C7
    net.source.c7 = net %>% dplyr::filter(source == 'F5_C7') %>% dplyr::filter(target %in% unique(epi@active.ident))
    #volcanoplot
    {
      volcano.dat = net.source.c7
      volcano.dat$significant = 'Stable'
      volcano.dat$significant[volcano.dat$ligand.logFC >=0.5] = 'up'
      volcano.dat$significant[volcano.dat$ligand.logFC <=-0.5] = 'down'

      volcano.dat1 = na.omit(volcano.dat)
      volcano.dat1 =volcano.dat1 %>% dplyr::distinct(interaction_name,.keep_all = T)
      volcano.dat1 =volcano.dat1 %>% dplyr::distinct(receptor,.keep_all = T)
      volcano.dat1 =volcano.dat1 %>% dplyr::distinct(ligand,.keep_all = T)
      
      volcano.dat1$interaction_name = as.character(volcano.dat1$interaction_name)
      
      pdf('/work/cwt/NACT_OV/3_cellchat/DEG_C7.pdf',width = 4.4,height = 3)
      ggplot(volcano.dat1,aes(ligand.logFC,-1*log10(ligand.pvalues),label = ifelse(abs(volcano.dat1$ligand.logFC)>0.75 ,interaction_name,'')))+
        geom_point(aes(color=significant))+theme_minimal()+
        geom_text_repel()+
        labs(title = "Differential LR Pairs of F5_C7",
             x = "Log Fold Change of Ligand",
             y = "-Log10(P-value)")+
        scale_color_manual(values = c("up" = "#D54062", "down" = "#2D5C7F", "Stable" = "grey"))
      dev.off()
    }
    #netplot
    {
      c7.up = subsetCommunication(cellchat, net = net.source.c7, datasets = "Post",ligand.logFC = 0.05, receptor.logFC = NULL)
      pair.c7.up = c7.up[,'interaction_name',drop = F]
      
      pdf('/work/cwt/NACT_OV/3_cellchat/F5_C7_postpre.pdf',width = 6.7,height = 5.5)
      netVisual_bubble(cellchat, sources.use = c(10), targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, 
                       title.name = "Increased signaling in Post-treatment", angle.x = 45, remove.isolate = F,
                       pairLR.use = pair.c7.up)
      dev.off()
    }
    
    p1 = VlnPlot(epi_fib,features = c('APP','CD74'),pt.size=0,stack = T)+scale_fill_manual(values = c('APP' = '#C1D8C3','CD74' = '#CD5C08'))
    p2 = VlnPlot(epi_fib,features = c('THBS1','CD47'),pt.size=0,stack = T)+scale_fill_manual(values = c('THBS1' = '#C8A1E0','CD47' = '#FF8C9E'))
    
    pdf('/work/cwt/NACT_OV/3_cellchat/vlnplot.pdf',width = 7,height = 5.6)
    p1+p2
    dev.off()
    
    
    #F8_CLEC3B
    {
      net.source.clec3b = net %>% dplyr::filter(source == 'F8_CLEC3B') %>% dplyr::filter(target %in% unique(epi@active.ident))
      
      clec3b.up = subsetCommunication(cellchat, net = net.source.clec3b, datasets = "Post",ligand.logFC = 0.05, receptor.logFC = NULL)
      pair.clec3b.up = clec3b.up[,'interaction_name',drop = F]
      
      netVisual_bubble(cellchat, pairLR.use = pair.clec3b.up, sources.use = c(6:13), targets.use = c(1:5), 
                       comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
      
      
      clec3b.down = subsetCommunication(cellchat, net = net.source.clec3b, datasets = "Pre",ligand.logFC = -0.05, receptor.logFC = NULL)
      pair.clec3b.down = clec3b.down[,'interaction_name',drop = F]
      netVisual_bubble(cellchat, pairLR.use = pair.clec3b.down, sources.use = c(13), targets.use = c(1:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
      
      
      #app-cd74 pair
      p1 = netVisual_bubble(cellchat, sources.use = c(13), targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, 
                            title.name = "Increased signaling in Post-treatment", angle.x = 45, remove.isolate = F,
                            pairLR.use = pair.clec3b.up)
      netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(13),  comparison = c(1, 2), max.dataset = 2, 
                       title.name = "Increased signaling in Post-treatment", angle.x = 45, remove.isolate = F)
      
      
      pdf('/work/cwt/NACT_OV/3_cellchat/app-cd74.pdf',width = 6.6,height = 1.8)
      p1
      dev.off()
      
      ##Volcanoplot of pair
      head(net.source.clec3b)
      
      volcano.dat = net.source.clec3b
      volcano.dat$significant = 'Stable'
      volcano.dat$significant[volcano.dat$ligand.logFC >=0.5] = 'up'
      volcano.dat$significant[volcano.dat$ligand.logFC <=-0.5] = 'down'
      
      table(volcano.dat$significant)
      volcano.dat1 = na.omit(volcano.dat)
      volcano.dat1 =volcano.dat1 %>% dplyr::distinct(interaction_name,.keep_all = T)
      
      getwd()
      volcano.dat1$interaction_name = as.character(volcano.dat1$interaction_name)
      
      pdf('/work/cwt/NACT_OV/3_cellchat/DEG_clec3b.pdf',width = 4.4,height = 3)
      ggplot(volcano.dat1,aes(ligand.logFC,-1*log10(ligand.pvalues),label = ifelse(volcano.dat1$significant != 'Stable',interaction_name,'')))+
        geom_point(aes(color=significant))+theme_minimal()+
        geom_text_repel()+
        labs(title = "Differential LR Pairs of F8_CLEC3B",
             x = "Log Fold Change of Ligand",
             y = "-Log10(P-value)")+
        scale_color_manual(values = c("up" = "#D54062", "down" = "#2D5C7F", "Stable" = "grey"))
      dev.off()
      
    }
    

     
     
     
     }
  
  
  
}