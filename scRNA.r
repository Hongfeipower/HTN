library(Seurat)
setwd("~/experiment/normal/RNA/data")
normal_rna = ReadMtx('UniqueAndMult-EM.mtx','barcodes.tsv','features.tsv')
setwd("~/experiment/HTN/RNA/data")
HTN_rna =ReadMtx('matrix.mtx','barcodes.tsv','features.tsv')
scRNAlist <- list()
scRNAlist[[1]]=CreateSeuratObject(normal_rna,project = 'normal',min.cells = 3,min.features = 200)
scRNAlist[[2]]=CreateSeuratObject(HTN_rna,project = 'HTN',min.cells = 3,min.features = 200)

scRNAlist[[1]][["percent.mt"]] = PercentageFeatureSet(scRNAlist[[1]], pattern = "^MT-")
scRNAlist[[2]][["percent.mt"]] = PercentageFeatureSet(scRNAlist[[2]], pattern = "^MT-")

VlnPlot(scRNAlist[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(scRNAlist[[2]], feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNAlist[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scRNAlist[[2]]=subset(scRNAlist[[2]],subset = nFeature_RNA>200 &nFeature_RNA<7500 & nCount_RNA<2700&percent.mt<2.4)
scRNAlist[[1]]=subset(scRNAlist[[1]],subset = nFeature_RNA>200 &nFeature_RNA<14000 & percent.mt<4.5)

scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

normal_top10 <- head(VariableFeatures(scRNAlist[[1]]), 10)
plot1 <- VariableFeaturePlot(scRNAlist[[1]])
plot2 <- LabelPoints(plot = plot1, points = normal_top10, repel = TRUE)
plot1 + plot2

HTN_top10 <- head(VariableFeatures(scRNAlist[[2]]), 10)
plot3 <- VariableFeaturePlot(scRNAlist[[2]])
plot4<- LabelPoints(plot = plot3, points = HTN_top10, repel = TRUE)
plot3+plot4

features <- SelectIntegrationFeatures(object.list = scRNAlist)

integration.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features)

integration.combined <- IntegrateData(anchorset = integration.anchors)

DefaultAssay(integration.combined) <- "integrated"
integration.combined <- ScaleData(integration.combined, verbose = FALSE)
integration.combined <- RunPCA(integration.combined, npcs = 30, verbose = FALSE)
integration.combined <- RunUMAP(integration.combined, reduction = "pca", dims = 1:30)
integration.combined <- FindNeighbors(integration.combined, reduction = "pca", dims = 1:30)
integration.combined <- FindClusters(integration.combined, resolution = 0.5)
pdf("integrate_RNA.pdf",width=16, height=8)
p1 <- DimPlot(integration.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(integration.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

pdf("integrate_RNA_split1.pdf",width=16, height=8)
DimPlot(integration.combined, reduction = "umap", split.by = "orig.ident",label = TRUE)
dev.off()
#DefaultAssay_integrated
all_markers <-FindAllMarkers(integration.combined,min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
markers_integrate=all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(markers_integrate,'top_markers_integrate.csv',row.names = FALSE)
write.csv(all_markers,'all_markers_intergrate.csv',row.names = FALSE)



#Cell annotation
library(SingleR)
load("~/experiment/integrate_RNA/hpca.se.RData")
#ref <- BlueprintEncodeData() 
clusters <- integration.combined@meta.data$seurat_clusters
testdata = GetAssayData(integration.combined, slot="data")
table(clusters)

cellpred <- SingleR(test = testdata,  
                    ref = hpca.se, 
                    labels = hpca.se$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels,
                      stringsAsFactors = F)





saveRDS(integration.combined,'integration.combined_no_celltype.rds')

# rename cell_type
Idents(integration.combined) <- "seurat_clusters"
integration.combined<-RenameIdents(integration.combined, `0` = "Hepatocytes type1", `1` = "Hepatocytes type2", `2` = "Endothelial",
                                 `3` = "Hepatocytes type3", `4` = "Hepatocytes type4", `5` = "Kupffer", `6` = "Stellate", `7` = "Hepatocytes type5", `8` = "Hepatocytes type6", `9` = "Hepatocytes type7",
                                 `10` = "T", `11` = "Hepatocytes type8", `12` = "Hepatocytes type9", `13` = "Hepatocytes type10", `14` = "Hepatocytes type11",
                              `15`="Hepatocytes type12",`16`='Stellate type2',`17`="B",`18`='Endothelial type2','19'='Erythroblast')

DimPlot(integration.combined, label = TRUE,split.by = "orig.ident")

saveRDS(integration.combined,'integration.combined_have_celltype.rds')


# Cell annotation based on manual annotation
library(Seurat)
integration.combined_have_celltype$new_type =integration.combined_have_celltype$celltype
Idents(integration.combined_have_celltype) <- "new_type"
integration<-RenameIdents(integration.combined_have_celltype, `Hepatocytes type1` = "Hepatocytes_1", `Hepatocytes type2` = "Hepatocytes_2", `Endothelial` = "Endothelial",
                                    `Hepatocytes type3` = "Hepatocytes_1", `Hepatocytes type4` = "Hepatocytes_2", `Kupffer` = "Kupffer", `Stellate` = "Stellate", `Hepatocytes type5` = "Hepatocytes_2", `Hepatocytes type6` = "Hepatocytes_1", `Hepatocytes type7` = "Cholangiocyte",
                                    `T` = "T_cell", `Hepatocytes type8` = "Hepatocytes_1", `Hepatocytes type9` = "Hepatocytes_1", `Hepatocytes type10` = "Hepatocytes_1", `Hepatocytes type11` = "Hepatocytes_2",
                                    `Hepatocytes type12`="Hepatocytes_1",`Stellate type2`='Stellate',`B`="Endothelial",`Endothelial type2`='Endothelial','Erythroblast'='Hepatocytes_2')
integration$new_type<-Idents(integration)
saveRDS(integration,'integration.rds')

all_markers <-FindAllMarkers(integration,min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
markers_integrate=all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(markers_integrate,'top_markers_newtype.csv',row.names = FALSE)
write.csv(all_markers,'all_markers_newtype.csv',row.names = FALSE)

integration$new_grop <- paste(Idents(integration), integration$orig.ident, sep = "_")
Idents(integration) <- "new_grop"
Hepatocytes_grup= FindMarkers(integration, ident.1 = 'Hepatocytes_2_HTN', ident.2 = 'Hepatocytes_2_normal', verbose = FALSE)
write.csv(Hepatocytes_grup,'Hepatocytes_2_grup_differential_gene.csv')
saveRDS(integration,'scRNA.rds')
