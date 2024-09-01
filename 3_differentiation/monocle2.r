{
  library(monocle)
  library(tidyverse)
  library(ggrastr)
  library(ggforce)
  
  getwd()
  setwd("/work/cwt/NACT_OV/6_differentiation")
    
  
  if(T){
    newimport <- function(otherCDS, import_all = FALSE) {
      if(class(otherCDS)[1] == 'Seurat') {
        requireNamespace("Seurat")
        data <- otherCDS@assays$RNA@counts
        
        if(class(data) == "data.frame") {
          data <- as(as.matrix(data), "sparseMatrix")
        }
        
        pd <- tryCatch( {
          pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
          pd
        }, 
        #warning = function(w) { },
        error = function(e) { 
          pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
          pd <- new("AnnotatedDataFrame", data = pData)
          
          message("This Seurat object doesn't provide any meta data");
          pd
        })
        
        # remove filtered cells from Seurat
        if(length(setdiff(colnames(data), rownames(pd))) > 0) {
          data <- data[, rownames(pd)]  
        }
        
        fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
        fd <- new("AnnotatedDataFrame", data = fData)
        lowerDetectionLimit <- 0
        
        if(all(data == floor(data))) {
          expressionFamily <- negbinomial.size()
        } else if(any(data < 0)){
          expressionFamily <- uninormal()
        } else {
          expressionFamily <- tobit()
        }
        
        valid_data <- data[, row.names(pd)]
        
        monocle_cds <- newCellDataSet(data,
                                      phenoData = pd, 
                                      featureData = fd,
                                      lowerDetectionLimit=lowerDetectionLimit,
                                      expressionFamily=expressionFamily)
        
        if(import_all) {
          if("Monocle" %in% names(otherCDS@misc)) {
            otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
            otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
            
            monocle_cds <- otherCDS@misc$Monocle
            mist_list <- otherCDS
            
          } else {
            # mist_list <- list(ident = ident) 
            mist_list <- otherCDS
          }
        } else {
          mist_list <- list()
        }
        
        if(1==1) {
          var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
          
        }
        monocle_cds@auxClusteringData$seurat <- mist_list
        
      } else if (class(otherCDS)[1] == 'SCESet') {
        requireNamespace("scater")
        
        message('Converting the exprs data in log scale back to original scale ...')    
        data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
        
        fd <- otherCDS@featureData
        pd <- otherCDS@phenoData
        experimentData = otherCDS@experimentData
        if("is.expr" %in% slotNames(otherCDS))
          lowerDetectionLimit <- otherCDS@is.expr
        else 
          lowerDetectionLimit <- 1
        
        if(all(data == floor(data))) {
          expressionFamily <- negbinomial.size()
        } else if(any(data < 0)){
          expressionFamily <- uninormal()
        } else {
          expressionFamily <- tobit()
        }
        
        if(import_all) {
          # mist_list <- list(iotherCDS@sc3,
          #                   otherCDS@reducedDimension)
          mist_list <- otherCDS 
          
        } else {
          mist_list <- list()
        }
        
        monocle_cds <- newCellDataSet(data,
                                      phenoData = pd, 
                                      featureData = fd,
                                      lowerDetectionLimit=lowerDetectionLimit,
                                      expressionFamily=expressionFamily)
        # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
        # monocle_cds@auxOrderingData$scran <- mist_list
        
        monocle_cds@auxOrderingData$scran <- mist_list
        
      } else {
        stop('the object type you want to export to is not supported yet')
      }
      
      return(monocle_cds)
    }
  }
  
  
  table(fib@active.ident)

  sname = subset(fib,downsample = 300)
```r

data <- as(as.matrix(sname@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sname@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())
```

sub <- newimport(sname)
sub <- estimateSizeFactors(sub) ### 计算Size_factor
sub <- estimateDispersions(sub) ### 计算Size_factor，评估离散度

head(pData(sub))
head(fData(sub))

diff_test_res <- differentialGeneTest(sub,fullModelFormulaStr = "~dcelltype")
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-100)) ###Inte the DE genes for construct paths
sub <- setOrderingFilter(sub, ordering_genes)
plot_ordering_genes(sub)
sub <- reduceDimension(
  sub,
  max_components = 2,
  method = 'DDRTree')

sub <- orderCells(sub,reverse = F) ##reverse= T?? 
saveRDS(sub,file = 'monocle_obj_fib.rds')
```

plot_cell_trajectory(sub, color_by = "dcelltype",show_branch_points = T)+
  scale_color_manual(values =fib.color)

plot_cell_trajectory(sub,color_by = 'Pseudotime', cell_size = 0.8,show_branch_points = T)



p2 = plot_cell_trajectory(sub, color_by = "dcelltype",show_branch_points = T)+
  scale_color_manual(values =fib.color)+NoLegend()+NoAxes()

p1 = plot_cell_trajectory(sub,color_by = 'Pseudotime', cell_size = 0.8,show_branch_points = T)+
  scale_colour_gradientn(colors = c('#CAE8D5','#3B6978','#204051'))+
  theme(legend.position = "bottom")+NoAxes()

p4 = plot_cell_trajectory(sub,color_by = 'treatment_phase', cell_size = 0.8,show_branch_points = F)+
  scale_color_manual(values = c('#2D5C7F','#984A59'))+NoLegend()+facet_wrap('~treatment_phase',nrow=1)+NoAxes();p4

p5= plot_cell_trajectory(sub,color_by = 'CA125_response', cell_size = 0.8,show_branch_points = F)+
  scale_color_manual(values = c('#2D5C7F','#984A59'))+NoLegend()+facet_wrap('~CA125_response',nrow=1)+NoAxes();p5

(p1+p2)/p4

ggsave(filename = 'monocle2_res_fib.pdf',width = 6,height = 5.8)


