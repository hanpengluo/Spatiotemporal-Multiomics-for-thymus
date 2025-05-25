# 00 package --------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(decontX)
library(DoubletFinder)
library(clustree)
library(tidyverse)
library(patchwork)
library(ggfun)

source('function_code.R')


# 01 combine --------------------------------------------------------------

thy_cite_all <- lapply(dir_list,FUN = function(project_name){
  thy_cite_rna <- readRDS(file.path(result_dir,file.path('02.CLEAN.RDS',paste0(project_name,'.clean.RDS'))))
  thy_cite_adt <- readRDS(file.path(result_dir,file.path('01.RAW.RDS',paste0(project_name,'.adt.raw.RDS'))))
  cell_match <- intersect(colnames(thy_cite_adt), colnames(thy_cite_rna))
  thy_cite_adt <- thy_cite_adt[, cell_match]
  seurat_cite <- thy_cite_rna[, cell_match]
  # Add ADT as an assay in Seurat object
  seurat_cite[['ADT']] <- CreateAssayObject(counts = thy_cite_adt)
  return(seurat_cite)
})
thy_cite_all <- merge(thy_cite_all[[1]],thy_cite_all[-1],add.cell.ids =dir_list)
thy_cite_all <- JoinLayers(thy_cite_all)


thy_cite_all <- seurat_reductions(thy_cite_all,res = 0.3)


# 02 CD19 clean -----------------------------------------------------------


DefaultAssay(thy_cite_all) <- 'ADT'
thy_cite_all_clean <- subset(thy_cite_all,CD19<0.5)


# 03 T cell clustered by CD4 and CD8a -------------------------------------

cite_adt_info <- FetchData(thy_cite_all_clean,vars = paste0('adt_',c("CD3","CD4","CD8a"))) %>%
  mutate(celltype = case_when(
    # DN
    adt_CD4 <= 0.35 & adt_CD8a <= 0.8 ~ 'DN',
    # DP
    adt_CD4 >= 0.4 & adt_CD8a >= 1 ~ 'DP',
    # CD4
    adt_CD4 >= 0.7 & adt_CD8a <= 0.8 ~ "CD4+T",
    # CD8
    adt_CD4 <= 0.3 & adt_CD8a >= 1.2 ~ "CD8+T",
    TRUE ~ "Other Cells"))
# add metadata with subtype of T
thy_cite_all_clean <- AddMetaData(thy_cite_all_clean,dplyr::select(cite_adt_info,celltype))
thy_cite_all_clean <- seurat_reductions(thy_cite_all_clean,res = 0.3)

thy_cite_all_clean_T <- subset(thy_cite_all_clean, celltype != 'Other Cells' ) %>% 
  seurat_reductions(res = 0.3)


# 04 cytotrace ------------------------------------------------------------

library(CytoTRACE2)
thy_cite_all_clean_T <- cytotrace2(input = thy_cite_all_clean_T,species = 'mouse',is_seurat = T,ncores = 12,seed = 123)
write.csv(thy_cite_all_clean_T@meta.data[,grep('CytoTRACE2',colnames(thy_cite_all_clean_T@meta.data),value = T)],
          'data_analysi/03.analysis/thy_cite_all_clean_T_cytotrace2.csv',row.names = T,quote = F)


# 05 sub population -------------------------------------------------------

# sub cell and wsnn
cite_cd4 <- seurat_cite_sub_f(thy_cite_all_clean_T,Celltype = "CD4+T",seq(0.01,0.1,0.01))
cite_cd8 <- seurat_cite_sub_f(thy_cite_all_clean_T,Celltype = "CD8+T",seq(0.01,0.1,0.01))
cite_dn <- seurat_cite_sub_f(thy_cite_all_clean_T,Celltype = "DN",seq(0.01,0.1,0.01))
cite_dp <- seurat_cite_sub_f(thy_cite_all_clean_T,Celltype = "DP",seq(0.01,0.1,0.01))
# create cell_number name
cite_cd8 <- iden_seurat(cite_cd8,ident_tag = 'wsnn_res.0.03')
cite_dp <- iden_seurat(cite_dp,ident_tag = 'wsnn_res.0.04')
cite_dn <- iden_seurat(cite_dn,ident_tag = 'wsnn_res.0.04')
cite_cd4 <- iden_seurat(cite_cd4,ident_tag = 'wsnn_res.0.05')


# 06 celltype annotation --------------------------------------------------

# DN -
Idents(sctype_dn) <- 'cell_number'
sctype_dn <- c("C_1"="DN3","C_2"="DN4","C_3"="DN1","C_4"="NKT","C_5"="Marcophage","C_6"="DC","C_7"="DN2","C_8"="B cell")
cite_dn@meta.data$sctype <- sctype_dn[as.character(cite_dn$cell_number)]
saveRDS(cite_dn,'data_analysi/03.analysis/cell.final.sctype.dn.RDS')

# DP -

Idents(sctype_dp) <- 'cell_number'
sctype_dp <- c("C_1"="DP2","C_2"="DP3","C_3"="DP stem","C_4"="DP1")
cite_dp@meta.data$sctype <- sctype_dp[as.character(cite_dp$cell_number)]
saveRDS(cite_dp,'data_analysi/03.analysis/cell.final.sctype.dp.RDS')

# CD4 -

Idents(cite_cd4) <- 'cell_number'
cite_cd4 <- FindSubCluster(cite_cd4,cluster = 'C_1',graph.name = 'wsnn',resolution = 0.2,algorithm = 3)
cite_cd4$cell_number <- cite_cd4$sub.cluster
cite_cd4$sub.cluster <- NULL
sctype_cd4 <- c("C_1_0"="Immature CD4 SP","C_1_1"="Mature CD4 SP","C_1_2"="Mature CD4 SP","C_1_3"="CD4 SP stem","C_2"="CD4 Treg")
cite_cd4@meta.data$sctype <- sctype_cd4[as.character(cite_cd4$cell_number)]
saveRDS(cite_cd4,'data_analysi/03.analysis/cell.final.sctype.cd4.RDS')

# CD8 -
Idents(cite_cd8) <- 'cell_number'
cite_cd8 <- FindSubCluster(cite_cd8,cluster = 'C_1',graph.name = 'wsnn',resolution = 0.2,algorithm = 3)
cite_cd8$cell_number <- cite_cd8$sub.cluster
cite_cd8$sub.cluster <- NULL
sctype_cd8 <- c("C_1_0"="Immature CD8 SP","C_1_1"="Mature CD8 SP","C_1_2"="Mature CD8 SP","C_2"="CD8 SP memory","C_3"="CD8 SP stem")
cite_cd8@meta.data$sctype <- sctype_cd8[as.character(cite_cd8$cell_number)]
saveRDS(cite_cd8,'data_analysi/03.analysis/cell.final.sctype.cd8.RDS')








