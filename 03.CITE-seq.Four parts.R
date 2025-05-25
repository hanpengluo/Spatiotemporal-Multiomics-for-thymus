# 00 package --------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)

source('function_code.R')


# 01 DN -------------------------------------------------------------------


# 01.1 import DN ---------------------------------------------------------------
# saveRDS(cite_dn,'data_analysi/03.analysis/cell.final.sctype.dn.RDS')
cite_dn <- readRDS('data_analysis/03.analysis/cell.final.sctype.dn.RDS')

color_dn_sctype <- c("DN3"="#8e0002","DN4"="#3F84AA","DN1"="#f704f6","NKT"="#ff9b94",
                     "Marcophage"="#c473a6","DC"="#17bd4f","DN2"="#f4b21f","B cell"="#999999",
                     "DN1_γδT"="#57C3F3","DN2_γδT"="#C5DEBA"
)
dn_result <- 'data_analysis/03.analysis/01.DN'

deg_f(cite_dn,ident_tag = 'sctype',deg_dir = file.path(dn_result,'DEG_dn8'),project_name = 'DN8')

# 01.2 γδT definition ---------------------------------------------------------------------

Idents(cite_dn)='sctype'
sctype_dn <- FetchData(cite_dn,vars = c('sctype', 
                                        grep('^Trac', rownames(cite_dn@assays$RNA), value = TRUE), 
                                        grep('^Trdc', rownames(cite_dn@assays$RNA), value = TRUE)), slot = 'data') %>%
  mutate(Trac_avg = rowMeans(dplyr::select(., starts_with("rna_Trac")), na.rm = TRUE),
         Trdc_avg = rowMeans(dplyr::select(., starts_with("rna_Trdc")), na.rm = TRUE),
         sctype = ifelse(Trac_avg == 0 & sctype=='DN1', 'DN1_γδT', 
                         ifelse(Trdc_avg > 0 & sctype=='DN2', 'DN2_γδT',sctype)))
cite_dn <- AddMetaData(cite_dn,dplyr::select(sctype_dn,sctype))

# 01.3 DEG for γδT ---------------------------------------------------------------------

deg_f(subset(cite_dn,sctype %in% c('DN1','DN1_γδT')),ident_tag = 'sctype',deg_dir = file.path(dn_result,'DEG_dn1'),project_name = 'DN1')
deg_f(subset(cite_dn,sctype %in% c('DN2', 'DN2_γδT')),ident_tag = 'sctype',deg_dir = file.path(dn_result,'DEG_dn2'),project_name = 'DN2')
deg_f(subset(cite_dn,sctype %in% c('DN1','DN1_γδT','DN2','DN2_γδT','DN3','DN4')),ident_tag = 'sctype',deg_dir = file.path(dn_result,'DEG_dn6'),project_name = 'DN6')

# 01.4 TF -------------------------------------------------------------------

deg_dn <- openxlsx::read.xlsx(file.path(dn_result,'DEG_dn6/Table.DN6.DEG.xlsx'))
cite_dn6 <- subset(cite_dn,sctype %in% c('DN1','DN1_γδT','DN2','DN2_γδT','DN3','DN4'))

tf_function(outdir = file.path(dn_result,'TF_dn6'),
            seurat_input = cite_dn6,
            genesubset = unique(deg_dn$gene),
            Org = 'mgi',myDatasetTitle = 'DN1',
            dbDir = '/mnt/d/01fasta/RcisTarget',
            #dbs = c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather','10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'),
            dbs = c('500bp'= 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                    '10kb' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'),
            core = 10,tf_method = 'GRNBoost',pyscenic_path='~/miniconda3/envs/pyscenic/bin/pyscenic',
            celltype_col = "sctype")

tf_result_function(outdir = file.path(dn_result,'TF_dn6'),
                   seurat_input = cite_dn6,
                   celltype_col='sctype')





# 02 DP -------------------------------------------------------------------



# 02.1 import DP-------------------------------------------------------------------

cite_dp <- readRDS('data_analysis/03.analysis/cell.final.sctype.dp.RDS')
dp_result <- 'data_analysis/03.analysis/02.DP'


# 02.2 DEG -------------------------------------------------------------------

deg_f(cite_dp,ident_tag = 'sctype',deg_dir = file.path(dp_result,'DEG_dp4'),project_name = 'DP4')
deg_f(subset(cite_dp,sctype %in% c('DP1','DP2','DP3')),ident_tag = 'sctype',deg_dir = file.path(dp_result,'DEG_dp3'),project_name = 'DP3')



# 02.3 TF -------------------------------------------------------------------

deg_dp <- openxlsx::read.xlsx(file.path(dp_result,'DEG_dp3/Table.DP3.DEG.xlsx'))
cite_dp3 <- subset(cite_dp,sctype %in% c('DP1','DP2','DP3'))
tf_function(outdir = file.path(dp_result,'TF_dp3'),
            seurat_input = cite_dp3,
            genesubset = unique(deg_dp$gene),
            Org = 'mgi',myDatasetTitle = 'DP1',
            dbDir = '/mnt/d/01fasta/RcisTarget',
            #dbs = c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather','10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'),
            dbs = c('500bp'= 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                    '10kb' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'),
            core = 10,tf_method = 'GRNBoost',pyscenic_path='~/miniconda3/envs/pyscenic/bin/pyscenic',
            celltype_col = "sctype")
tf_result_function(outdir = file.path(dp_result,'TF_dp3'),
                   seurat_input = cite_dp3,
                   celltype_col='sctype')

# 02.4 TF activity ------------------------------------------------------------

scenicOptions <- readRDS("int/1.5_scenicOptions.Rds")
cite_dp3 <- subset(cite_dp, sctype %in% c('DP1','DP2','DP3'))
cellInfo <- cite_dp3@meta.data[, c('orig.ident','nCount_RNA','nFeature_RNA','sctype')]
colnames(cellInfo)[4] <- 'CellType'
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
rownames(regulonAUC) <- gsub("_extended", "", rownames(regulonAUC))
cite_dp3[['regulonAUC']] <- CreateAssayObject(counts = getAUC(regulonAUC))

set.seed(123)
subset_seurat <- subset(cite_dp3, cells = unlist(lapply(
  split(colnames(cite_dp3), cite_dp3$sctype), 
  function(x) sample(x, min(500, length(x)))
)))

regulonAUC_subset <- as.data.frame(subset_seurat@assays$regulonAUC$counts)
rownames(regulonAUC_subset) <- gsub("\\(.*?\\)", "(+)", rownames(regulonAUC_subset))
ht_auc_scale <- t(scale(t(as.matrix(regulonAUC_subset))))

cell_order <- order(factor(subset_seurat$sctype, levels = c("DP1", "DP2", "DP3")))
ht_auc_scale_ordered <- ht_auc_scale[, cell_order]

tf_groups <- list(
  DP1 = c("Foxp4 (+)"),
  DP2 = c("Nfatc3 (+)", "Nfkb1 (+)", "Nfkb2 (+)", "Rel (+)", "Relb (+)", "Junb (+)", "Fos (+)","Rorc (+)"),
  DP3 = c("Nfatc1 (+)", "Runx3 (+)", "Trp53 (+)", "Stat6 (+)")
)

group_colors <- c("DP1" = "#FF1493", "DP2" = "#8E44AD", "DP3" = "#2E8B57")
tf_colors <- unlist(lapply(names(tf_groups), function(x) rep(group_colors[x], length(tf_groups[[x]]))))
names(tf_colors) <- unlist(tf_groups)

row_anno <- rowAnnotation(
  mark = anno_mark(
    at = which(rownames(ht_auc_scale_ordered) %in% names(tf_colors)),
    labels = rownames(ht_auc_scale_ordered)[rownames(ht_auc_scale_ordered) %in% names(tf_colors)],
    labels_gp = gpar(
      col = tf_colors[rownames(ht_auc_scale_ordered)[rownames(ht_auc_scale_ordered) %in% names(tf_colors)]],
      fontface = "bold",
      fontsize = 10
    ),
    link_width = unit(5, "mm"),
    padding = unit(1, "mm")
  ),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  matrix = ht_auc_scale_ordered,
  top_annotation = col_anno,
  name = "Regulon activity (Z-score)",
  col = colorRamp2(c(-4, 0, 4), c("white", "#FCF8DE", "#253177")),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  right_annotation = row_anno,
  use_raster = TRUE,
  row_dend_side = "left",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8),
    legend_width = unit(5, "cm"),
    legend_height = unit(0.1, "cm")
  )
)




# 03 CD4 ------------------------------------------------------------------



# 03.1 CD4 import ----------------------------------------------------------------------
cite_cd4 <- readRDS('data_analysis/03.analysis/cell.final.sctype.cd4.RDS')
cd4_result <- 'data_analysis/03.analysis/03.CD4'

# 03.2 CD4 DEG ----------------------------------------------------------------------

# CD4 all sctype
deg_f(seurat_obj = cite_cd4,ident_tag = 'sctype',deg_dir = file.path(cd4_result,'DEG_cd4_4'),project_name = 'CD4')

# CD4 Immature vs Mature
cite_cd4sub <- subset(cite_cd4,sctype %in% c("Immature CD4 SP", "Mature CD4 SP"))
deg_f(cite_cd4sub,ident_tag = 'sctype',deg_dir = file.path(cd4_result,'DEG_cd4_mature_immature'),project_name = 'CD4_mature_immature')

# CD4 Treg vs cd4T
cite_cd4@meta.data$treg <- ifelse(cite_cd4$sctype=='CD4 Treg','CD4 Treg','CD4 T')
deg_f(seurat_obj = cite_cd4,ident_tag = 'treg',deg_dir = file.path(cd4_result,'DEG_Treg'),project_name = 'Treg')


# 03.3 TF for mature & immature -------------------------------------------

deg_cd4_mature_immature <- openxlsx::read.xlsx(file.path(cd4_result,'DEG_cd4_mature_immature/Table.CD4_mature_immature.DEG.xlsx'))

cite_cd4mature <- subset(cite_cd4,sctype %in% c('Immature CD4 SP','Mature CD4 SP'))
tf_function(outdir = file.path(cd4_result,'TF_cd4_mature_immature'),
            seurat_input = cite_cd4mature,
            genesubset = unique(deg_cd4_mature_immature$gene),
            Org = 'mgi',myDatasetTitle = 'CD4',
            dbDir = '/mnt/d/01fasta/RcisTarget',
            #dbs = c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather','10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'),
            dbs = c('500bp'= 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                    '10kb' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'),
            core = 10,tf_method = 'GRNBoost',pyscenic_path='~/miniconda3/envs/pyscenic/bin/pyscenic',
            celltype_col = "sctype")
tf_result_function(outdir = file.path(cd4_result,'TF_cd4_mature_immature'),
                   seurat_input = cite_cd4mature,
                   celltype_col='sctype')



# 03.4 DEG for DP3 ------------------------------------------------

thy_cite_dp3 <- subset(thy_cite_all_clean_T,sctype =='DP3')
exp_data <- FetchData(thy_cite_dp3, vars = c("Gata3", "Runx3"))
DefaultAssay(thy_cite_dp3) <- 'RNA'
exp_data <- FetchData(thy_cite_dp3,vars = c("Gata3", "Runx3")) %>% 
  mutate(grtype = case_when(
    # gata3+runx-
    Gata3 > 0 & Runx3 == 0 ~ 'gata3+',
    # runx3
    Gata3 == 0 & Runx3 > 0 ~ 'runx3+',
    TRUE ~ "Other Cells"))
table(exp_data$grtype)
thy_cite_dp3 <- AddMetaData(thy_cite_dp3,dplyr::select(exp_data,grtype))

deg_f(seurat_obj = subset(thy_cite_dp3 ,subset=grtype %in% c('gata3+','runx3+')),ident_tag = 'grtype',deg_dir = file.path(cd4_result,'Gata3_Runx3'),project_name = 'Gata3_Runx3')



# 03.5 validation --------------------------------------------------------------------
# publication data GSE186078
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186078
# quality control for each sample same as 01.QC.R
pub_dir <- 'Thymus_pub/GSE186078_RAW'
seurat_thymus <- readRDS(file.path(pub_dir,'02.thymus.reduction.RDS'))
info <- openxlsx::read.xlsx('Thymus_pub/GSE186078_RAW/sample_info.xlsx')
sample_detail <- paste0(info$Genotype,'-',info$Location,'-',info$FACS_sort)
names(sample_detail) <- info$gsm
seurat_thymus@meta.data$sample_detail <- sample_detail[as.character(seurat_thymus$orig.ident)]

# 1. CD4+CD8+ T cell definition -
cite_adt_info <- FetchData(seurat_thymus,vars = paste0('adt_',c("CD3","CD4","CD8a"))) %>% 
  mutate(celltype = case_when(
    # DN
    adt_CD4 <= 3 & adt_CD8a <= 3 ~ 'DN',
    # DP
    adt_CD4 >= 3.5 & adt_CD8a >= 3.5 ~ 'DP',
    # CD4
    adt_CD4 >= 4 & adt_CD8a <= 3 ~ "CD4+T",
    # CD8
    adt_CD4 <= 3 & adt_CD8a >= 4 ~ "CD8+T",
    TRUE ~ "Other Cells"))
# 2. add celltype to meta data
seurat_thymus <- AddMetaData(seurat_thymus,dplyr::select(cite_adt_info,celltype))

# 3. mature process CD24 and CD62L
p_24_62_gene_heat_function <- function(rna_name){
  p <- FetchData(subset(seurat_thymus,subset = celltype =='CD4+T' & sample_detail %in% c('B6-Berkeley-CD5+TCRB+','MHCII-/--Berkeley-CD5+TCRB+','OT-I-BioLegend-NA')),
                 vars = unique(c('sample_detail','CD24','CD62L',rna_name)),layer = 'data') %>% 
    ggplot(aes(x=adt_CD62L,y=adt_CD24))+
    geom_point(aes_string(color=rna_name,size=rna_name),alpha=1)+
    scale_size(range = c(0.5,1))+
    geom_density_2d(show.legend = F,alpha=1,color='black',bins = 12)+
    labs(x='Surface protein CD62L',y='Surface protein CD24',title=rna_name)+
    xlim(c(0,6))+scale_x_continuous(breaks = 1:6)+ylim(c(0,5))+
    facet_wrap(~sample_detail,
               labeller = labeller(sample_detail = c(
                 'B6-Berkeley-CD5+TCRB+' = 'B6',
                 'MHCII-/--Berkeley-CD5+TCRB+' = 'MHCII-/-',
                 'OT-I-BioLegend-NA' = 'OT-I')))+
    scale_color_gradientn(colors = c("blue", "yellow", "red"))+
    cowplot::theme_cowplot(font_size = 18)+
    theme(legend.title = element_blank(),
          axis.line = element_blank(),strip.background = element_rect(fill = "white", color = "white"),
          strip.text = element_text(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
    )
  return(p)
}
p_rna_S1pr1 <- p_24_62_gene_heat_function('rna_S1pr1')
p_rna_S1pr1
ggsave(p_rna_S1pr1,filename = file.path(pub_dir,'Figure/Figure.cd62l_cd24_rna_S1pr1.pdf'),width = 10,height = 4.5)


# 4. CD4+CD8+ T cell surface protein
p_thy_single_function <- function(sample_id,color_name,title_name){
  p <- FetchData(seurat_thymus,vars = unique(c('sample','sample_detail',paste0('adt_',c("CD3","CD4","CD8a")),'rna_Gata3','rna_Runx3'))) %>% 
    dplyr::filter(sample_detail==sample_id) %>% 
    ggplot(aes(adt_CD4,adt_CD8a))+
    geom_point(color=color_name,size=0.15,alpha=1)+
    geom_density_2d(show.legend = F,alpha=1,color='black',bins = 8)+
    annotate("rect", xmin = 4, ymax = 3, xmax = max(cite_adt_info$adt_CD4), ymin = 0, color = 'black', fill = NA, linewidth = 1,linetype=2) +
    labs(title=title_name,color='',x="Surface protein: CD4",y="Surface protein: CD8a")+
    scale_color_gradientn(colors = c("blue", "yellow", "red"))+
    theme(plot.title = element_text(hjust = 0.5))+
    cowplot::theme_cowplot(font_size = 18)+
    theme(legend.position = 'top',
          legend.title = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3))
  return(p)
}
p_b6_thy <- p_thy_single_function('B6-Berkeley-CD5+TCRB+',color_all[1],'B6')
p_b6_thy
p_mhc1_thy <- p_thy_single_function('MHCI-/--Berkeley-CD5+TCRB+',color_all[2],'MHCI-/-')
p_mhc1_thy
p_ot1_thy <- p_thy_single_function('OT-I-BioLegend-NA',color_all[3],'OT-I')
p_ot1_thy
p_mhc2_thy <- p_thy_single_function('MHCII-/--Berkeley-CD5+TCRB+',color_all[4],'MHCII-/-')
p_mhc2_thy
p_ot2_thy <- p_thy_single_function('OT-II-BioLegend-NA',color_all[5],'OT-II')
p_ot2_thy

p <- (p_b6_thy | p_mhc1_thy | p_mhc2_thy )/
  (p_ot1_thy | p_ot2_thy | plot_spacer())
p
ggsave(p,file = file.path(pub_dir,'Figure/Figure.thymus.CD4&8.pdf'),width = 10,height = 7)






# 04 CD8 ------------------------------------------------------------------



# 04.1 CD8 import ---------------------------------------------------------

cite_cd8 <- readRDS('data_analysis/03.analysis/cell.final.sctype.cd8.RDS')
cd8_result <- 'data_analysis/03.analysis/04.CD8'

# 04.2 CD8 DEG ------------------------------------------------------------

# deg all
deg_f(seurat_obj = cite_cd8,ident_tag = 'sctype',deg_dir = file.path(cd8_result,'DEG_cd8_4'),project_name = 'CD8')

# DEG cd8_mature_immature
cite_cd8sub <- subset(cite_cd8,sctype %in% c("Immature CD8 SP", "Mature CD8 SP"))
deg_f(cite_cd8sub,ident_tag = 'sctype',deg_dir = file.path(cd8_result,'DEG_CD8_mature_immature'),project_name = 'CD8_mature_immature')

# deg stem
cite_cd8@meta.data$stem <- ifelse(cite_cd8$sctype=='CD8 SP stem','CD8 SP stem','CD8T')
deg_f(seurat_obj = cite_cd8,ident_tag = 'stem',deg_dir = file.path(cd8_result,'DEG_cd8_stem'),project_name = 'CD8_sp_stem')

# 04.3 TF for mature & immature -------------------------------------------
deg_cd8_mature_immature <- openxlsx::read.xlsx(file.path(cd8_result,'DEG_CD8_mature_immature/Table.CD8_mature_immature.DEG.xlsx'))
cite_cd8mature <- subset(cite_cd8,sctype %in% c('Immature CD8 SP','Mature CD8 SP'))
tf_function(outdir = file.path(cd8_result,'TF_cd8_mature_immature'),
            seurat_input = cite_cd8mature,
            genesubset = unique(deg_cd8_mature_immature$gene),
            Org = 'mgi',myDatasetTitle = 'CD8',
            dbDir = '/mnt/d/01fasta/RcisTarget',
            #dbs = c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather','10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'),
            dbs = c('500bp'= 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                    '10kb' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'),
            core = 10,tf_method = 'GRNBoost',pyscenic_path='~/miniconda3/envs/pyscenic/bin/pyscenic',
            celltype_col = "sctype")
tf_result_function(outdir = file.path(cd8_result,'TF_cd8_mature_immature'),
                   seurat_input = cite_cd8mature,
                   celltype_col='sctype')

# 04.4 TF for stem ------------------------------------------------------------
deg_cd8_stem <- openxlsx::read.xlsx(file.path(cd8_result,'DEG_cd8_stem/Table.CD8_sp_stem.DEG.xlsx'))
cite_cd8stem <- subset(cite_cd8,sctype %in% c('CD8 SP stem','CD8T'))
tf_function(outdir = file.path(cd8_result,'TF_cd8_stem'),
            seurat_input = cite_cd8stem,
            genesubset = unique(deg_cd8_stem$gene),
            Org = 'mgi',myDatasetTitle = 'CD8',
            dbDir = '/mnt/d/01fasta/RcisTarget',
            #dbs = c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather','10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'),
            dbs = c('500bp'= 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                    '10kb' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'),
            core = 10,tf_method = 'GRNBoost',pyscenic_path='~/miniconda3/envs/pyscenic/bin/pyscenic',
            celltype_col = "sctype")
tf_result_function(outdir = file.path(cd8_result,'TF_cd8_stem'),
                   seurat_input = cite_cd8stem,
                   celltype_col='sctype')
# 04.5 stem gene enrichment ------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
deg_stem <- openxlsx::read.xlsx(file.path(cd8_result,'DEG_cd8_stem/Table.CD8_sp_stem.DEG.xlsx'))
deg_gene <- bitr(deg_stem$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Mm.eg.db')
kegg <- enrichKEGG(gene = deg_gene$ENTREZID,organism = 'mmu',pvalueCutoff = 0.05)
kegg <- DOSE::setReadable(kegg,OrgDb = 'org.Mm.eg.db',keyType = 'ENTREZID')
kegg_result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_result$Description)
top10 <- kegg_result %>%
  arrange(p.adjust) %>%  
  head(10) %>%
  mutate(Description = factor(Description, levels = rev(Description)))
p_kegg_stem <- ggplot(top10, aes(x = Count, y = Description)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),fill = "grey98", alpha = 0.3, color = NA) +
  geom_segment(aes(x = 0, xend = Count, y = Description, yend = Description),color = "grey70", linewidth = 0.6) +
  geom_point(aes(size = Count, color = -log10(p.adjust)),alpha = 0.9, shape = 19) +
  viridis::scale_color_viridis(option = "D", end = 0.9, name = "-log10(p.adj)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count",breaks = scales::pretty_breaks(n = 4)) +
  labs(x = "Gene Count",y = NULL,title = "Top 10 KEGG Pathways for Stem Cell") +
  cowplot::theme_cowplot(font_size = 10)
ggsave(p_kegg_stem,filename = file.path(cd8_result,'Figure/Figure.kegg_stem.pdf'),width = 5,height = 4)


# 04.6 cell cycle pathway ------------------------------------------------
library(ggh4x)
library(ggfun)

genes <- unlist(str_split(kegg_result$geneID[kegg_result$Description=='Cell cycle'], '/'))
gene_label_df <- data.frame(
  gene = c(grep('Cc', genes, value = TRUE),
           grep('Cdk', genes, value = TRUE),
           grep('Ch', genes, value = TRUE),
           str_subset(unlist(str_split(kegg_result$geneID[1], '/')), "Trp|My")),
  lable = c(rep("Cc", length(grep('Cc', genes, value = TRUE))),
            rep("Cdk", length(grep('Cdk', genes, value = TRUE))),
            rep("Ch", length(grep('Ch', genes, value = TRUE))),
            rep("other", length(str_subset(unlist(str_split(kegg_result$geneID[1], '/')), "Trp|My"))))
)
gene_label_df <- merge(gene_label_df,deg_stem,by='gene') %>% dplyr::filter(!duplicated(gene))
plot_cycle_pathway <- ggplot(gene_label_df)+
  geom_bar(aes(x=-log10(p_val),y=interaction(gene,lable),fill=lable),stat="identity",show.legend = F)+
  geom_text(aes(x=1,y=interaction(gene,lable),label=gene),hjust="left",color="#000000",show.legend = F)+
  geom_point(aes(x=130,y=interaction(gene,lable),size=avg_log2FC,fill=lable),shape=21)+
  geom_text(aes(x=130,y=interaction(gene,lable),label=round(avg_log2FC,0)),show.legend = F)+
  scale_x_continuous(expand=expansion(mult=c(0,0.2)))+
  scale_fill_manual(values=c("#74c476","#41b6c4","#f46d43","#9e9ac8"))+
  scale_size(range=c(3,6),guide=guide_legend(override.aes=list(fill="#000000")))+
  guides(y="axis_nested",fill=guide_legend(reverse=T))+labs(x="-log10(pvalue)",y="Cell cycle gene class",size='LogFC')+
  theme_bw()+theme(
    ggh4x.axis.nestline.y=element_line(size=3,color=c("#74c476","#41b6c4","#f46d43","#9e9ac8")),
    ggh4x.axis.nesttext.y=element_text(colour=c("#74c476","#41b6c4","#f46d43","#9e9ac8")),
    panel.border=element_rect(linewidth=1),
    axis.text=element_text(color="#000000",size=14),
    axis.text.y=element_text(color=rep(c("#41ae76","#225ea8","#fc4e2a","#88419d"),each=10)),
    axis.text.y.left=element_blank(),
    axis.ticks.length.y.left=unit(10,"pt"),
    axis.ticks.y.left=element_line(color=NA),
    axis.title=element_text(color="#000000",size=15),
    plot.title=element_text(color="#000000",size=20,hjust=0.5)
  )+guides(fill=FALSE)
ggsave(plot_cycle_pathway,filename = file.path(cd8_result,'Figure/Figure.kegg_stem_cell_cycle.pdf'),width = 5,height = 4)


# 04.7 DP SP --------------------------------------------------------------

thy_cite_dpsp <- subset(thy_cite_all_clean_T,sctype %in% c('DP1','DP2','DP3','Immature CD4 SP','Mature CD4 SP','Immature CD8 SP','Mature CD8 SP'))
sce <- thy_cite_dpsp
DefaultAssay(sce) <- "RNA"
sce[["UMAP"]] <- sce[["wnn.umap"]]
sce <- SeuratWrappers::as.cell_data_set(x = sce,assay = 'RNA')
sce <- monocle3::cluster_cells(cds = sce,reduction_method = 'UMAP')
sce <- monocle3::learn_graph(sce, use_partition = F )
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                        (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
sce <- monocle3::order_cells(sce, root_pr_nodes=myselect(sce,select.classify = 'sctype',my_select = "DP1"))

pseudotime_cds <- monocle3::pseudotime(sce,reduction_method = 'UMAP') %>% data.frame(pseudotime=.)
thy_cite_dpsp <- AddMetaData(thy_cite_dpsp,pseudotime_cds)

plot_gene_expression <- function(df, gene, exclude_cd4 = c('Immature CD4 SP', 'Mature CD4 SP'), exclude_cd8 = c('Immature CD8 SP', 'Mature CD8 SP'), color_all_sctype) {
  # DP->CD4SP
  p1 <- df %>% 
    dplyr::filter(! sctype %in% exclude_cd8) %>% 
    ggplot(aes(x = pseudotime, y = !!sym(gene))) +
    geom_point(aes(fill = sctype), shape = 21, alpha = 1) +
    geom_smooth(color = 'red', size = 2) +
    labs(x = "Pseudotime", y = paste0(gene, " Expression"), fill = "DP->CD4SP") +
    cowplot::theme_cowplot(font_size = 8) +
    scale_fill_manual(values = color_all_sctype) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  # DP->CD8SP
  p2 <- df %>% 
    dplyr::filter(! sctype %in% exclude_cd4) %>% 
    ggplot(aes(x = pseudotime, y = !!sym(gene))) +
    geom_point(aes(fill = sctype), shape = 21, alpha = 1) +
    geom_smooth(color = 'red', size = 2) +
    labs(x = "Pseudotime", y = paste0(gene, " Expression"), fill = "DP->CD8SP") +
    cowplot::theme_cowplot(font_size = 8) +
    scale_fill_manual(values = color_all_sctype) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  # 返回两个图的组合
  p <- (p1 | p2)
  return(p)
}

p_runx3 <- plot_gene_expression(df, gene = "rna_Runx3", color_all_sctype = color_all_sctype)
p_gata3 <- plot_gene_expression(df, gene = "rna_Gata3", color_all_sctype = color_all_sctype)
p_zbtb7b <- plot_gene_expression(df, gene = "rna_Zbtb7b", color_all_sctype = color_all_sctype)
p_smad7 <- plot_gene_expression(df, gene = "rna_Smad7", color_all_sctype = color_all_sctype)

p_monocle3_dpsp <- p_runx3 / p_gata3 / p_smad7 / p_zbtb7b
p_monocle3_dpsp
ggsave(p_monocle3_dpsp,filename = file.path(cd8_result,'Figure/Figure.monocle3_dpsp.pdf'),width = 5,height = 7)
