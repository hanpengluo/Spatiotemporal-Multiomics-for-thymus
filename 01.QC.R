
# 00 package --------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(decontX)
library(DoubletFinder)
library(tidyverse)
library(patchwork)


# 01. quality conyrol ---------------------------------------------------------------------

dir_list <- gsub(pattern = '_gex','',list.files('03.count',recursive = F,pattern = '_gex$'))

for(i in dir_list){
  cellranger_out <- file.path(paste0('03.count','/',i),'outs/filtered_feature_bc_matrix')
  thy_cite <- Read10X(cellranger_out)
  thy_cite_rna <- thy_cite$`Gene Expression`
  thy_cite_adt <- thy_cite$`Antibody Capture`
  adt_annotation <- read.csv('Totalseq_B_Mouse_tag_cellranger.csv')
  rownames(thy_cite_adt) <- adt_annotation$tag[match(adt_annotation$name,rownames(thy_cite_adt))]
  
  project_name=i
  mtname="^mt-"
  result_dir <- 'data_analysi'
  if(!dir.exists(result_dir)) dir.create(result_dir)
  
  sce <- SingleCellExperiment(list(counts=thy_cite_rna))
  sce <- decontX(sce)
  plot_ambient_rna <- plotDecontXContamination(sce)
  df_con_rate <- data.frame(decontX=sce$decontX_contamination)
  rownames(df_con_rate) <- colnames(sce)
  seurat_obj <- CreateSeuratObject(counts = round(decontXcounts(sce)), assay = 'RNA',project = project_name, min.cells = 3, min.features = 200)
  seurat_obj <- AddMetaData(seurat_obj,df_con_rate)
  seurat_obj@meta.data$log10GenesPerUMI <- log10(seurat_obj@meta.data$nFeature_RNA)/log10(seurat_obj@meta.data$nCount_RNA)
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = mtname, col.name = "percent.mt")
  # threshold
  nFeature_RNA_min=200
  nFeature_RNA_max=6000
  con_rate=0.9
  GenesPerUMI=0.8
  seuratdf <- subset(seurat_obj,
                     subset = nFeature_RNA > nFeature_RNA_min & 
                       nFeature_RNA < nFeature_RNA_max & 
                       decontX < con_rate & 
                       log10GenesPerUMI > GenesPerUMI)
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","decontX","log10GenesPerUMI"), cols = 'blue3',ncol = 5,
                group.by = 'orig.ident',pt.size = 0) & labs(x = 'RAW') + theme(text=element_text(size=6))
  p2 <- VlnPlot(seuratdf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","decontX","log10GenesPerUMI"),cols = 'red3', ncol = 5,
                group.by = 'orig.ident',pt.size = 0) & labs(x = 'QC')+ theme(text=element_text(size=6))
  p1/p2
  rds_path <- file.path(result_dir,'01.RAW.RDS')
  if (!dir.exists(rds_path)) {
    dir.create(rds_path)
  }
  saveRDS(seurat_obj,file = file.path(rds_path,paste0(project_name,'.rna.raw.RDS')))
  saveRDS(thy_cite_adt,file = file.path(rds_path,paste0(project_name,'.adt.raw.RDS')))
  # workflow
  seuratdf <- seuratdf %>% 
    NormalizeData(verbose = FALSE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
    RunPCA(npcs = 30, verbose = FALSE) 
  
  #Select number of pca
  pct <- seuratdf [['pca']]@stdev / sum(seuratdf [['pca']]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  
  seuratdf <-  seuratdf %>%  RunUMAP(reduction = "pca", dims = 1:pcs)%>% 
    FindNeighbors(reduction = "pca", dims = 1:pcs) %>% 
    FindClusters(resolution = 0.5)
  
  # doublelet
  sweep.res.list_kidney <- paramSweep(seuratdf, PCs = 1:pcs,num.cores = 1)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  ## PK
  mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  print(paste0('Best parameter of PK is ',mpK))
  #Homotypic Doublet Proportion Estimate
  annotations <- seuratdf@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(seuratdf)*8*1e-6
  print(paste0('The Doublet Proportion is ',DoubletRate))
  nExp_poi <- round(DoubletRate*nrow(seuratdf@meta.data))
  # doublet
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  pn=0.25
  # singlet
  seuratdf <- doubletFinder(seuratdf, PCs = 1:pcs, pN = pn, pK = mpK, nExp = nExp_poi, reuse.pANN = F, sct = F)
  seuratdf <- doubletFinder(seuratdf, PCs = 1:pcs, pN = pn, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  seuratdf[[paste('pANN',pn,mpK,nExp_poi,sep = '_')]] <- NULL
  seuratdf[[paste('pANN',pn,mpK,nExp_poi.adj,sep = '_')]] <- NULL
  names(seuratdf@meta.data)[which(names(seuratdf@meta.data) %in% c(paste('DF.classifications',pn,mpK,nExp_poi,sep = '_'),paste('DF.classifications',pn,mpK,nExp_poi.adj,sep = '_'))
  )] <- c('dfpoi','dfpoi2')
  seuratdf@meta.data$doublefinder <- ifelse(seuratdf@meta.data$dfpoi=='Doublet' & seuratdf@meta.data$dfpoi2=='Doublet','Doublet-High Confidience',
                                            ifelse(seuratdf@meta.data$dfpoi=='Singlet' & seuratdf@meta.data$dfpoi2=='Singlet','Singlet','Doublet-Low Confidience'))
  # table(seuratdf@meta.data$doublefinder)
  p3 <- DimPlot(seuratdf, reduction = "umap", group.by ="doublefinder",cols =c("black","red","gold"),pt.size = 1.3)
  seuratdf_clean <- subset(seuratdf,doublefinder=='Singlet')
  # save Seurat
  rds_clean_path <- file.path(result_dir,'02.CLEAN.RDS')
  if (!dir.exists(rds_clean_path)) {
    dir.create(rds_clean_path)
  }
  saveRDS(seuratdf,file = file.path(rds_clean_path,paste0(project_name,'.Singlet.RDS')))
  raw_counts <- GetAssayData(seuratdf_clean, assay = "RNA", layer = "counts")
  new_seurat_obj <- CreateSeuratObject(counts = raw_counts, project = project_name)
  saveRDS(new_seurat_obj,file = file.path(rds_clean_path,paste0(project_name,'.clean.RDS')))
  # report
  qc_data <- data.frame(
    Title = project_name,
    Dataset = c('raw_data', 'qc_raw_data' ,'qc_feature_data', 'qc_singlet_data'),
    Ngenes = c(dim(thy_cite_rna)[1],dim(seurat_obj)[1],dim(seuratdf)[1], dim(new_seurat_obj)[1]),
    Ncells = c(dim(thy_cite_rna)[2], dim(seurat_obj)[2],dim(seuratdf)[2], dim(new_seurat_obj)[2]),
    QCgene = c(0, dim(thy_cite_rna)[1] - dim(seurat_obj)[1],
               dim(seurat_obj)[1] - dim(seuratdf)[1], 
               dim(seuratdf)[1] - dim(new_seurat_obj)[1]),
    QCcell = c(0, dim(thy_cite_rna)[2] - dim(seurat_obj)[2],
               dim(seurat_obj)[2] - dim(seuratdf)[2],
               dim(seuratdf)[2] - dim(new_seurat_obj)[2]),
    QCcell_precentage=c(0,paste0(100*round((dim(thy_cite_rna)[2] - dim(seurat_obj)[2])/dim(thy_cite_rna)[2],4),'%'),
                        paste0(100*round((dim(seurat_obj)[2] - dim(seuratdf)[2])/dim(thy_cite_rna)[2],4),'%'),
                        paste0(100*round((dim(seuratdf)[2] - dim(new_seurat_obj)[2])/dim(thy_cite_rna)[2],4),'%'))
  )
  write.csv(qc_data, file = file.path(rds_clean_path, paste0(project_name, '.qc.csv')), row.names = FALSE, quote = FALSE)
  p <- ((p1/p2) |( plot_ambient_rna/p3))+patchwork::plot_layout(widths = c(3,1))
  ggsave(file = file.path(rds_clean_path, paste0(project_name, '.qc.png')), plot = p, width = 18, height = 6)
  print(paste0(project_name,' is Done'))
}



