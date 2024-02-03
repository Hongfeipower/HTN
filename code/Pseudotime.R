library(monocle)
data <- as(as.matrix(scRNA@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~ new_grop")


ordering_genes <- row.names (subset(diff_test_res, qval < 0.1)) ## 

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 3,
                        num_dim = 20,
                        method = 'DDRTree') # DDRTree·½Ê½HSMM <- orderCells(HSMM)
HSMM <- orderCells(HSMM)


colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


a1 <- plot_cell_trajectory(HSMM, color_by = "State") + scale_color_manual(values = colour)

a2 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime") 

a3<- plot_cell_trajectory(HSMM, color_by = "new_type") + scale_color_manual(values = colour)

a4 <- plot_cell_trajectory(HSMM, color_by = "new_grop")  + scale_color_manual(values = colour)

a<-HSMM@reducedDimS

a=t(a)

b<-HSMM@phenoData@data$Pseudotime

b=cbind(rownames(a),b)
write.csv(b,'track_cell_Pseudotime.csv')
