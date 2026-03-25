# Convert rds to h5ad file ####
library(Seurat)
library(SeuratDisk)

epi <- readRDS('R_objects/THCA_scrna_epithelial.rds')
epi$annotation_epi_scanpy <- as.character(epi_sub3$annotation_epi)
SaveH5Seurat(epi, filename = "palantir/THCA_scrna_epithelial.h5Seurat")
Convert("palantir/THCA_scrna_epithelial.h5Seurat", dest = "h5ad")

#initial and terminal cell decision using monocle2####
library(monocle)
library(SoupX)
library(MAST)
epi <- readRDS('R_objects/THCA_scrna_epithelial.rds')
Idents(epi) <- 'annotation_epi'
epi_mnn_DEG_Top300 <- quickMarkers(epi@assays$RNA@counts, epi$annotation_epi, N=300)
DefaultAssay(epi) <- 'RNA'
cds <- as.CellDataSet(epi)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- monocle::setOrderingFilter(cds, unique(epi_mnn_DEG_Top300$gene))
cds <- monocle::reduceDimension(cds, method = 'DDRTree')
cds <- monocle::orderCells(cds)
cds$annotation_epi <- factor(cds$annotation_epi, levels = c('Follicular cell', 'PTC type1', 'PTC type2', 'ATC'))
row.names(epi@meta.data)[which.min(cds$Pseudotime)] # PTC3_AGCGTATCAATGAATG-1
row.names(epi@meta.data)[which.max(cds$Pseudotime)] # ATC3_CTTAGGACAGTATGCT-1
saveRDS(cds, "monocle/THCA_scrna_epithelial_monocle2_cds.rds")


# TF expression trend along the trajectory ####
interest_regulons_significant_clustered <- read.csv('scenic_epi_sub3_macbook/exploring_data_including_follicular/thca_by_cancertype_clustering_average/20240417_thca_by_cancertype_interesting_clustering_wilcox_mcquitty_regulons_km.csv',
                                                    row.names = 1)
tf <- as.character(do.call(rbind.data.frame, strsplit(interest_regulons_significant_clustered$x, " "))[[1]])
tf <- gsub('_extended','',tf)
write.csv(tf, 'palantir/THCA_scrna_epithelial_scenic_significant_tf_including_follicular.csv')


## Fig. 4D - module score####
tf_cluster <- read.csv('palantir/20240429_thca_scrna_epi_sub3_palantir_TF_including_follicular_gene_trend_cluster_km.csv', row.names = 1)
gene_cluster1 <- tf_cluster %>% filter(Cluster == '1') %>% pull(Genes)
gene_cluster4 <- tf_cluster %>% filter(Cluster == '4') %>% pull(Genes)
gene_cluster0 <- tf_cluster %>% filter(Cluster == '0') %>% pull(Genes)
gene_cluster2 <- tf_cluster %>% filter(Cluster == '2') %>% pull(Genes)
gene_cluster3 <- tf_cluster %>% filter(Cluster == '3') %>% pull(Genes)

epi_sub3 <- AddModuleScore(epi_sub3, features = list(gene_cluster1), name="TF_cluster1")
epi_sub3 <- AddModuleScore(epi_sub3, features = list(gene_cluster4), name="TF_cluster4")
epi_sub3 <- AddModuleScore(epi_sub3, features = list(gene_cluster0), name="TF_cluster0")
epi_sub3 <- AddModuleScore(epi_sub3, features = list(gene_cluster2), name="TF_cluster2")
epi_sub3 <- AddModuleScore(epi_sub3, features = list(gene_cluster3), name="TF_cluster3")

tf_ms <- epi_sub3@meta.data %>% dplyr::select(TF_cluster21, TF_cluster11, TF_cluster01, TF_cluster31, TF_cluster41)

write.csv(tf_ms, 'palantir/epi_sub3_TF_including_follicular_cluster_module_score.csv', row.names = T)


## Fig. 4E - Pathways enriched in each of the TF clusters####
library(enrichR)
regulons <- readRDS('scenic_epi_sub3_macbook/int/2.6_regulons_asGeneSet.Rds')

tf_cluster <- read.csv('palantir/20240429_thca_scrna_epi_sub3_palantir_TF_including_follicular_gene_trend_cluster_km.csv', row.names = 1)
gene_cluster0 <- tf_cluster %>% filter(Cluster == '0') %>% pull(Genes)
gene_cluster1 <- tf_cluster %>% filter(Cluster == '1') %>% pull(Genes)
gene_cluster2 <- tf_cluster %>% filter(Cluster == '2') %>% pull(Genes)
gene_cluster3 <- tf_cluster %>% filter(Cluster == '3') %>% pull(Genes)
gene_cluster4 <- tf_cluster %>% filter(Cluster == '4') %>% pull(Genes)

regulons_unlist <- unlist(regulons)
#C0
cluster0_regulon <- c()
for (i in gene_cluster0) {
  genesets <- regulons_unlist[grep(paste0(i,'[1-9]'),names(regulons_unlist))]
  cluster0_regulon <- c(cluster0_regulon, genesets)
}
cluster0_regulon <- unique(cluster0_regulon)
length(cluster0_regulon)

#C1
cluster1_regulon <- c()
for (i in gene_cluster1) {
  genesets <- regulons_unlist[grep(paste0(i,'[1-9]'),names(regulons_unlist))]
  cluster1_regulon <- c(cluster1_regulon, genesets)
}
cluster1_regulon <- unique(cluster1_regulon)

#C2
cluster2_regulon <- c()
for (i in gene_cluster2) {
  genesets <- regulons_unlist[grep(paste0(i,'[1-9]'),names(regulons_unlist))]
  cluster2_regulon <- c(cluster2_regulon, genesets)
}
cluster2_regulon <- unique(cluster2_regulon)

#C3
cluster3_regulon <- c()
for (i in gene_cluster3) {
  genesets <- regulons_unlist[grep(paste0(i,'[1-9]'),names(regulons_unlist))]
  cluster3_regulon <- c(cluster3_regulon, genesets)
}
cluster3_regulon <- unique(cluster3_regulon)

#C4
cluster4_regulon <- c()
for (i in gene_cluster4) {
  genesets <- regulons_unlist[grep(paste0(i,'[1-9]'),names(regulons_unlist))]
  cluster4_regulon <- c(cluster4_regulon, genesets)
}
cluster4_regulon <- unique(cluster4_regulon)



dbs <- c("MSigDB_Hallmark_2020")
clusters <- c(1,4,0,2,3)
df_enriched_flt_all_regulons <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Term","Adjusted.P.value", "Combined.Score", "Genes", "Cluster"))
for(i in clusters) {
  geneset <- get(paste0('cluster',i,'_regulon'))
  df_enriched <- as.data.frame(enrichr(geneset, dbs))
  colnames(df_enriched) <- gsub("MSigDB_Hallmark_2020.", "", colnames(df_enriched))
  df_enriched$Genes <- gsub(";", ", ", df_enriched$Genes)
  df_enriched_flt <- df_enriched %>% select(Term, Adjusted.P.value, Combined.Score, Genes) %>% filter(Adjusted.P.value < 0.05) %>% arrange(desc(Combined.Score))
  df_enriched_flt$Cluster <- i
  df_enriched_flt_all_regulons <- rbind(df_enriched_flt_all_regulons,df_enriched_flt)
}

library(openxlsx)
write.xlsx(df_enriched_flt_all_regulons, "palantir/20240429_thca_scrna_epi_sub3_palantir_TF_including_follicular_gene_trend_cluster_regulons_erichr_km.xlsx")


# Fig. 4G - correlation between potency score and CEBPB expression level ####
library(CytoTRACE2)
library(Seurat)
library(dplyr)
library(openxlsx)
library(RColorBrewer)
library(viridis)

epi <- readRDS('R_objects/THCA_scrna_epithelial.rds')
cytotrace2_result <- cytotrace2(epi, species = 'human', is_seurat = TRUE, slot_type = 'counts', batch_size = 10000, parallelize_models = TRUE, ncores = 4, seed = 1234)
cyto_res <- cytotrace2_result@meta.data %>% select(CytoTRACE2_Score,CytoTRACE2_Potency,CytoTRACE2_Relative,preKNN_CytoTRACE2_Score,preKNN_CytoTRACE2_Potency)

write.xlsx(cyto_res, 'cytotrace/20251010_thca_scrna_cytotrace_potency_km.xlsx', row.names = T)
cyto_res <- read.xlsx('cytotrace/20251010_thca_scrna_cytotrace_potency_km.xlsx', rowNames = T)
epi <- AddMetaData(epi,cyto_res)

cebpb_potency <- data.frame(
  cebpb = epi@assays$RNA@data['CEBPB',],
  potency = epi$CytoTRACE2_Score,
  group = epi$annotation_epi)

cebpb_potency$group <- factor(cebpb_potency$group, levels = c('Follicular cell', 'PTC type1', 'PTC type2', 'ATC'))

ggplot(cebpb_potency, aes(x = cebpb, y = potency)) +
  geom_point(aes(shape = group, color = group)) +
  scale_color_manual(values = c('#EA8331', '#7CAE00', '#35A2FF', '#C00000')) +
  geom_smooth(method = 'lm', color ='red') +
  annotate("text", x = 0.5, y = 0.6, col = "black", size = 10,
           label = paste("r = ", signif(cor(cebpb_potency$cebpb, cebpb_potency$potency),3)))+
  theme_classic() +
  labs(x = 'CEBPB expression level', y = 'Potency score') +
  theme(axis.title.x = element_text(size = 27, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 27, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 27, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 27, family = 'Arial', colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 27, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank()) &
  guides(shape=guide_legend(override.aes = list(size = 7)))


