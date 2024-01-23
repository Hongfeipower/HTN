setwd("/experiment/normal/ATAC/data")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 128) 
inputFiles='fragments.tsv.gz'
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = 'normal',
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

setwd("/experiment/integrate_ATAC")
ArrowFiles=c('HTN.arrow','normol.arrow')

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)

# dimensionality reduction and clustering

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- saveArchRProject(ArchRProj = proj)


# match scRNA and scATAC
setwd("~/experiment/integrate_ATAC")
library(ArchR)
library(Seurat)
set.seed(1)
addArchRThreads(threads = 128)
proj= loadArchRProject('HemeTutorial')

#find scRNA HTN and normal cellname
Idents(scRNA)<-'orig.ident'

rnaHTN<- WhichCells(scRNA, idents = 'HTN')
rnanormal<-WhichCells(scRNA, idents = 'normal')

# alignement of background informatiaon for sample types
groupList <- SimpleList(
  HTN = SimpleList(
    ATAC = proj$cellNames[proj$Sample=='HTN'],
    RNA = rnaHTN
  ),
  normal = SimpleList(
    ATAC = proj$cellNames[proj$Sample=='normol'],
    RNA = rnanormal
  )
)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = scRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "new_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)


pal <- paletteDiscrete(values=unique(scRNA$new_type))

p2 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "predictedGroup", 
  pal = pal
)

markerGenes<- c("G6PC",'PCK1',"PECAM1",'LDB2',"MARCO", "ACTA2", "KRT7", 'CD2')

projHeme3 <- addImputeWeights(proj)

p1 <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3)
)


p2 <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3)
)


p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 4), p2c))


plotPDF(plotList = p2, 
        name = "3-1_marke_GeneScoreMatrix.pdf", 
        ArchRProj = projHeme3, 
        addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
#CM = as.matrix(confusionMatrix(projHeme5$Clusters, projHeme5$predictedGroup))
#write.table(CM,'CM.csv')



labelOld <- rownames(cM)

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew, oldLabels = labelOld)
projHeme3$Cluster_type=paste0(projHeme3$Sample,'_',projHeme3$Clusters2)

p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2")
p2 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Cluster_type")
plotPDF(p1, name = "4-Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2, name = "5-Plot-UMAP-Remap-Clusters-type.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3", load = FALSE)
projHeme3= loadArchRProject('Save-ProjHeme3')

hpe1_normal = projHeme3$cellNames[projHeme3$Cluster_type=='normol_Hepatocytes_1']
hpe1_HTN = projHeme3$cellNames[projHeme3$Cluster_type=='HTN_Hepatocytes_1']

hpe2_normal= projHeme3$cellNames[projHeme3$Cluster_type=='normol_Hepatocytes_2'] 
hpe2_HTN = projHeme3$cellNames[projHeme3$Cluster_type=='HTN_Hepatocytes_2']

#step 4 call peak by MACS2
projHeme3 = loadArchRProject('Save-ProjHeme3')


markerGenes=c("PCCA", "CYP3A4", "PPARGC1A", "UPB1","SDS")
p <- plotBrowserTrack(
  ArchRProj = projHeme3, 
  groupBy = "Cluster_type",
  useGroups =c('HTN_Hepatocytes_1','normol_Hepatocytes_1'),
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 2000
)

grid::grid.newpage()
grid::grid.draw(p$PCCA)

plotPDF(plotList = p, 
        name = "6-Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projHeme3, 
        addDOC = FALSE, width = 5, height = 5)



addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- projHeme3@genomeAnnotation$genome 
BSgenome <- eval(parse(text = genome)) 
BSgenome <- validBSgenome(BSgenome) 

projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Cluster_type")
pathToMacs2 <- findMacs2()

projHeme4 <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Cluster_type", 
  pathToMacs2 = pathToMacs2
)
peak_information =getPeakSet(projHeme4)
write.csv(peak_information,'peak_information.csv')

# step 5 add peak
projHeme5 <- addPeakMatrix(projHeme4)
getAvailableMatrices(projHeme5)

# get Marker Peaks

markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Cluster_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# find up peak in one group

HTN_Hep <- markerPlot(seMarker = markersPeaks, name = "HTN_Hepatocytes_1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
normal_Hep <-markerPlot(seMarker = markersPeaks, name = "normal_Hepatocytes_1", cutOff = "FDR <= 0.1 & Log2FC >= 0.5", plotAs = "MA")



# setp 6get peak and peak background
a =getPeakSet(projHeme5)
b= getBgdPeaks(projHeme5)

write.table(a,sep='\t','all_tyep_peak.bed')
write.csv(b@rowRanges,'peak_bgd.csv')

# step 7 find marker peak

markerList=getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)
p <- plotBrowserTrack(
  ArchRProj = projHeme5, 
  groupBy = "Cluster_type", 
  geneSymbol = markerGenes,
  features =markerList[c('HTN_Hepatocytes_1','normol_Hepatocytes_1')],
  upstream = 50000,
  downstream = 2000,
  useGroups =c('HTN_Hepatocytes_1','normol_Hepatocytes_1')
)

grid::grid.newpage()
grid::grid.draw(p$SDS)
plotPDF(p, name = "6-2-Plot-Tracks-With-Features", width = 5, height = 8, ArchRProj = projHeme5, addDOC = FALSE)

saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-ProjHeme5", load = FALSE)


markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "HTN_Hepatocytes_1"
)

#step 8 add motif 
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

