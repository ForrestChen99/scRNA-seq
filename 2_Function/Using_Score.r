process_diff_combined <- function(Diff_HM, dcelltype, selected_hallmark) {
  # Set names for the list
  names(Diff_HM) <- dcelltype
  
  # Combine the list into a single dataframe
  Diff_combined <- dplyr::bind_rows(Diff_HM)
  rownames(Diff_combined) <- NULL
  
  # Modify the Pathway column
  Diff_combined <- Diff_combined %>%
    mutate(Pathway = gsub("HALLMARK-", "", terms))
  
  # Modify the Pathway names
  Diff_combined$Pathway <- gsub('-', '_', Diff_combined$Pathway)
  
  # Filter and factor the Diff_combined dataframe
  Diff_combined_selected <- dplyr::filter(Diff_combined, Pathway %in% selected_hallmark$Name)
  Diff_combined_selected$Pathway <- factor(Diff_combined_selected$Pathway, levels = selected_hallmark$Name)
  
  return(Diff_combined_selected)
}
  


  library(msigdbr)
  library(fgsea)
  library(GSVA)
  HM = gmtPathways('/work/cwt/NACT_OV/1_clustering/genelist/h.all.v2023.2.Hs.symbols.gmt')
  
    
  #selected_hallmark = readxl::read_xlsx(path = '/work/cwt/NACT_OV/1_clustering/genelist/HALLMARK_CATEGORY.xlsx')
  #selected_hallmark = arrange(selected_hallmark,Category)
  
  source('/work/cwt/ID8/2_function/scores.R')##YH DU
  HMscore <- score_cells(seur=epi, names=HM, combine_genes='mean', 
                         groups=NULL, group_stat='mean', cells.use=NULL)
  
  HMmatrix = as.matrix(HMscore) %>% t
  HMmatrix = apply(HMmatrix, 2, function(x)signif(x,digits = 3))
  colnames(HMmatrix) = rownames(epi@meta.data)
  
  HMseuobj = CreateSeuratObject(counts = HMmatrix,data = HMmatrix,meta.data = epi@meta.data)
  
  epi@assays$HM = HMseuobj@assays$RNA
  
 
  epi$Metastatic = ifelse(epi$position == 'primary','Primary','Metastatic')
  table(epi$Metastatic)
  
  dcelltype = unique(epi$dcelltype)
  
  
  DefaultAssay(epi) = 'HM'
  
  Diff_HM_treatment = lapply(dcelltype,function(i){
    sub = subset(epi,subset = dcelltype == i)
    DEfeatures =FindMarkers(sub,ident.1 = 'post',ident.2 = 'naive',group.by = 'treatment_phase',
                            min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
    DEfeatures$terms = rownames(DEfeatures)
    DEfeatures$celltype = i
    return(DEfeatures)
    print(i)
  })
  
  Diff_HM_metastatic = lapply(dcelltype,function(i){
    sub = subset(epi,subset = dcelltype == i)
    DEfeatures =FindMarkers(sub,ident.1 = 'Metastatic',ident.2 = 'Primary',group.by = 'Metastatic',
                            min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
    DEfeatures$terms = rownames(DEfeatures)
    DEfeatures$celltype = i
    return(DEfeatures)
    print(i)
  })
  
  Diff_HM_RESPONSE = lapply(dcelltype,function(i){
    sub = subset(epi,subset = dcelltype == i)
    DEfeatures =FindMarkers(sub,ident.1 = 'CR',ident.2 = 'PR',group.by = 'CA125_response',
                            min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
    DEfeatures$terms = rownames(DEfeatures)
    DEfeatures$celltype = i
    return(DEfeatures)
    print(i)
  })
  

  
  plot.data.treatment = process_diff_combined(Diff_HM_treatment,dcelltype,selected_hallmark)
  plot.data.metastatic = process_diff_combined(Diff_HM_metastatic,dcelltype,selected_hallmark)
  plot.data.response = process_diff_combined(Diff_HM_RESPONSE,dcelltype,selected_hallmark)
  
  p1 = ggplot(plot.data.treatment,aes(y = celltype,x = Pathway))+
    geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj+10^-100)))+
    scale_color_gradient2(low = "#769FCD", mid = "white", high = "#FC5185", midpoint = 0, 
                          limits = c(-1, 1))+
    labs(title = 'HALLMARK:Post vs Naive')+
    scale_size_continuous(range=c(1,6),breaks=seq(0,30,10),name="FDR")+
    theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
          panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=45,vjust = 0.5,hjust = 0.5,color="black"),
          axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+coord_flip()+
    guides(size = guide_legend(nrow = 2));p1
  
  p2 = ggplot(plot.data.metastatic,aes(y = celltype,x = Pathway))+
    geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj+10^-100)))+
    scale_color_gradient2(low = "#769FCD", mid = "white", high = "#FC5185", midpoint = 0, 
                          limits = c(-0.75, 0.75))+
    labs(title = 'HALLMARK:Metastatic vs Primary')+
    scale_size_continuous(range=c(1,6),breaks=seq(0,30,10),name="FDR")+
    theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
          panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=45,vjust = 0.5,hjust = 0.5,color="black"),
          axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+coord_flip()+
    guides(size = guide_legend(nrow = 2));p2
  
  p3 = ggplot(plot.data.response,aes(y = celltype,x = Pathway))+
    geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj+10^-100)))+
    scale_color_gradient2(low = "#769FCD", mid = "white", high = "#FC5185", midpoint = 0, 
                          limits = c(-0.5,0.5))+
    labs(title = 'HALLMARK:CR vs NR')+
    scale_size_continuous(range=c(1,6),breaks=seq(0,30,10),name="FDR")+
    theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
          panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=45,vjust = 0.5,hjust = 0.5,color="black"),
          axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+coord_flip()+
    guides(size = guide_legend(nrow = 2));p3
  
  cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol = 3)
  getwd()
  ggsave(filename = '/work/cwt/NACT_OV/2_function/hallmark.epi.pdf',width = 18.3,height = 13.2)
