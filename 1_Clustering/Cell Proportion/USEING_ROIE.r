##ROIE
{
  source('/work/cwt/ID8/1_clustering/ROIE.R')
  meta = tumor_ident@meta.data
  summary = table(meta[,c('dcelltype','group')])
  roe = as.data.frame(ROIE(summary))
  roe$cluster = rownames(roe)
  rownames(roe) = NULL
  res = data.frame()
  res = rbind(res,roe)
  summary(roe)

  heatmapdata = gather(roe,group,value,-c(cluster))
  str(heatmapdata)
  table(meta$bigclass)
  heatmapdata$group = factor(heatmapdata$group,levels = c('ID8-wt','ID8-PPK1'))
  heatmapdata$cluster = factor(heatmapdata$cluster,levels = c('TC1','TC2','TC3','TC4',
                                                              'TC5'))
  ggplot(heatmapdata,aes(x = group,y=cluster,fill = value))+
    geom_tile()+scale_fill_gradient(high = "#D7301F",low = "#FFF7EC")+
    geom_tile()+
    geom_text(aes(label = round(value,2)))+theme_transparent()+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5,face = "bold"))+ylab("")+xlab("")+
    theme(axis.text.y = element_text(face = "bold"))
  ggsave(path = "/work/cwt/ID8/1_clustering/",filename = "tumor_cell_compositio.pdf",width = 3,height = 5)
}