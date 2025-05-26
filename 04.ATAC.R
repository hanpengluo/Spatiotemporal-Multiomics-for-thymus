# 00 package --------------------------------------------------------------

library(Seurat)
library(JASPAR2020)
library(Signac)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(cicero)
library(SeuratWrappers)
library(patchwork)

future::plan("multisession", workers = 10) 
options(future.globals.maxSize = 100 * 1024 ^ 3)
source('function_code.R')


# 01 data import ----------------------------------------------------------------------

primary_input(metadata_file = 'ATAC/outs/singlecell.csv',
              atacpath='ATAC/outs/filtered_peak_bc_matrix.h5',
              fragpath="ATAC/outs/fragments.tsv.gz",
              outname="ATAC/01.thymus_atac.RDS")

# 02 quality control -------------------------------------------------------------

cite_thy_atac <- readRDS('ATAC/01.thymus_atac.RDS')
cite_thy_atac <- NucleosomeSignal(object = cite_thy_atac)
cite_thy_atac$nucleosome_group <- ifelse(cite_thy_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = cite_thy_atac, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
cite_thy_atac <- TSSEnrichment(cite_thy_atac)
cite_thy_atac <- subset(
  x = cite_thy_atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
cite_thy_atac <- RunTFIDF(cite_thy_atac)
cite_thy_atac <- FindTopFeatures(cite_thy_atac, min.cutoff = 'q0')
cite_thy_atac <- RunSVD(object = cite_thy_atac)
cite_thy_atac <- RunUMAP(object = cite_thy_atac,reduction = 'lsi',dims = 2:30)
cite_thy_atac <- FindNeighbors(object = cite_thy_atac,reduction = 'lsi',dims = 2:30)
cite_thy_atac <- FindClusters(object = cite_thy_atac,algorithm = 3,resolution = 1,verbose = FALSE)


# 03 gene activation score ----------------------------------------------------------

gene.activities <- GeneActivity(cite_thy_atac)
# add the gene activity matrix to the Seurat object as a new assay
cite_thy_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
cite_thy_atac <- NormalizeData(object = cite_thy_atac,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(cite_thy_atac$nCount_RNA))


# 04 mapquery -------------------------------------------------------------

thy_cite_all_clean_T <- FindVariableFeatures(object = thy_cite_all_clean_T,nfeatures = 5000)

transfer.anchors2 <- FindTransferAnchors(
  reference = thy_cite_all_clean_T,
  query = cite_thy_atac,
  reduction = 'cca', # cca rpca
  dims = 1:30
)

cite_thy_atac <- TransferData(
  anchorset = transfer.anchors1,
  reference = thy_cite_all_clean_T,
  query = cite_thy_atac,
  refdata = list(
    sctype = "sctype",
    sctype2 = "sctype2",
    predicted_ADT = "ADT"),
  weight.reduction = cite_thy_atac1[['lsi']],
  dims = 2:30
)

# flow plot

p_predictedADT <- FeatureScatter(cite_thy_atac,feature1 = 'CD4',feature2 = 'CD8a',group.by = 'predicted.sctype',cols = color_all_sctype)$data %>%
  ggplot(aes(predictedADT_CD4,predictedADT_CD8a))+
  scattermore::geom_scattermore(aes(color=colors),pointsize=4,alpha=1,pixels=c(1200,1200),interpolate=TRUE)+
  geom_density_2d(show.legend = F,alpha=1,color='black',bins = 10)+
  scale_color_manual(values = color_all_sctype)+
  labs(color='')+
  cowplot::theme_cowplot()
ggsave(p_predictedADT,filename = 'ATAC/Figures/predictedADT.pdf',width = 6,height = 6)


# 05 call peak ------------------------------------------------------------

peaks_all <- CallPeaks(
  object =  cite_thy_atac,
  assay = 'peaks',group.by = 'predicted.sctype2',
  macs2.path='~/miniconda3/envs/scell/bin/macs2'
)

peaks_all <- keepStandardChromosomes(peaks_all, pruning.mode = 'coarse')
peaks_all <- subsetByOverlaps(x = peaks_all,ranges = blacklist_mm10_unified, invert = TRUE)
frags <- CreateFragmentObject(
  path = 'ATAC/outs/fragments.tsv.gz',
  cells = colnames(cite_thy_atac),
  genome = BSgenome.Mmusculus.UCSC.mm10)

macs_count <- FeatureMatrix(fragments = frags,
                            features = peaks_all,
                            cells = colnames(cite_thy_atac))

cite_thy_atac[['ATAC']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'ATAC/outs/fragments.tsv.gz',
  annotation = annotations
)




# 06 accessibility ----------------------------------------------------------------------

cds <- as.cell_data_set(x = cite_thy_atac,assay = 'ATAC')
cds <- make_cicero_cds(cds,reduced_coordinates = reducedDims(cds)$UMAP)
# get the chromosome sizes from the Seurat object
genome <- seqlengths(cite_thy_atac)
# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome) %>% dplyr::filter(chr %in% paste0('chr',c(1:19,'X','Y','M')))
# run cicero
conns <- run_cicero(cds, genomic_coords = genome.df, sample_num = 100)
# Find cis-co-accessible networks (CCANs)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(cite_thy_atac) <- links


# 07 motif ----------------------------------------------------------------

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
cite_thy_atac <- AddMotifs(
  object = cite_thy_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
cite_thy_atac <- RunChromVAR(object = cite_thy_atac,genome = BSgenome.Mmusculus.UCSC.mm10)


# 08 LinkPeaks ----------------------------------------------------------------------

cite_thy_atac <-  Signac::RegionStats(object = cite_thy_atac,genome = BSgenome.Mmusculus.UCSC.mm10,assay = 'ATAC')
cite_thy_atac <- LinkPeaks(cite_thy_atac,peak.assay = 'ATAC',expression.assay = 'RNA',pvalue_cutoff = 0.05,genes.use = c('Rag1','Rag2'),distance = 5e+05)

link_atac_rag12 <- as.data.frame(Links(cite_thy_atac2@assays$ATAC))

ranges.show <- StringToGRanges('chr2-101552279-101552890')
ranges.show$color <- c("orange")

p <- CoveragePlot(
  object = cite_thy_atac2,
  # region = c("Rag1","Rag2"),
  region = 'chr2-101540000-101660000',
  peaks = F,
  group.by = 'predicted.sctype',
  region.highlight = ranges.show,
  idents = c('DN1','DN2','DN3','DN4','DP1','DP2','DP3'),
  extend.upstream = 1000,
  extend.downstream = 0
) & 
  theme(text = element_text(size=12),strip.placement = 'outside') & 
  scale_fill_manual(values = color_all_sctype)
p
ggsave(p,filename = 'ATAC/Figures/Figure.Rag1_Rag2_coverage.pdf',width = 8,height = 4)

