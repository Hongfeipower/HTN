library('tidyverse')
library(scRNAtoolVis)

data = read.csv('all_markers_intergrate.csv')
# name=c("Hepatocytes type1", "Hepatocytes type2",  "Endothelial",
#       "Hepatocytes type3", "Hepatocytes type4",  "Kupffer",  "Stellate",  "Hepatocytes type5",  "Hepatocytes type6", "Hepatocytes type7",
#        "T",  "Hepatocytes type8",  "Hepatocytes type9",  "Hepatocytes type10",  "Hepatocytes type11",
#        "Hepatocytes type12",'Stellate type2',"B",'Endothelial type2','Erythroblast')

for (i in c(0:19))
{
  data[data$cluster==i,"cluster"]=paste0('C',i)
}

#Hepatocytes= paste0('Hepatocytes type',c(1:12))

Hepatocytes= paste0('C',c(0,1,3,4,7,8,9,c(11:15)))
Hepatocytes_data = data[data$cluster %in% Hepatocytes, ]
`%!in%` <- Negate(`%in%`)
other_data= data[data$cluster %!in% Hepatocytes, ]


p1=jjVolcano(Hepatocytes_data, tile.col = corrplot::COL2('PuOr', 15)[3:15], size  = 2,topGeneN = 10,cluster.order=Hepatocytes)
p2=jjVolcano(other_data, tile.col = corrplot::COL2('PuOr', 15)[3:15], size  = 2,topGeneN = 10,cluster.order=c('C2','C18','C6','C16','C5','C10','C17','C19'))


Hepatocytes_gene_set=list()


for(i in Hepatocytes)
{
    cluster_data=Hepatocytes_data[Hepatocytes_data$cluster==i,]
    cluster_data_index =  order(cluster_data$avg_log2FC, decreasing = TRUE)[1:20]
    top_gene=cluster_data[ cluster_data_index,'gene']
    Hepatocytes_gene_set[[i]]=top_gene
} 


calculateJaccardSimilarity <- function(set1, set2) {
  # 计算交集的大小
  intersection_size <- length(intersect(set1, set2))
  
  # 计算并集的大小
  union_size <- length(union(set1, set2))
  
  # 计算Jaccard相似性
  jaccard_similarity <- intersection_size / union_size
  
  # 返回Jaccard相似性
  return(jaccard_similarity)
}


hot_data = diag(12)
for( i in c(2:12))
{
  print(c('i',i)) 
  for( j in c(1:(i-1)))
  {
    #print(c(i,j))
    hot_data[i,j]=calculateJaccardSimilarity(Hepatocytes_gene_set[[i]],Hepatocytes_gene_set[[j]])
    hot_data[j,i]=hot_data[i,j]
    
   
  }
}

rownames(hot_data)<-Hepatocytes
colnames(hot_data)<-Hepatocytes
pheatmap(hot_data, 
         cluster_rows = FALSE, 
         cluster_cols = TRUE,  
         main = "Clustered Heatmap", 
         fontsize = 8,  
         cellwidth = 15, 
         cellheight = 15,  
          display_numbers = TRUE
)

##find marker
Hepatocyte_marker=c("ALB", "AFP", "AHSG", "TTR", "SERPINA1", "APOA1", "GLUL", "ASL", "CYP2E1", "ASS1", "PCK1", "G6PC", "FABP1", "RBP4", "FGG", "FGA", "FGB", "APOC3", "APOA4", "GSTA1")
Endothelial_marker =c("PECAM1", "STAB2", "KDR", "OIT3", "IL1A", "F8", "GNG11", "LDB2", "PPAP2B", "GPR116", "FCN3", "HES1", "CRHBP", "STAB1", "DNASE1L3", "DLC1", "CLEC4G")
Cholangiocyte_marker=c('KRT8','KRT19','KRT7','SPP1','OPN3','sox9','HNF1','BPROM1')
Kupffer_marker=c('FCGR3','CD68','MARCO','IRF7','SPIC','CD14','TLR4',
                 'FCGR1A','FCGR1G','MSR1','CETP','CSF1R','TLF9','CLEC4F')

Stellate_marker=c('COL1A1','DCN','ACTA2','GFAP','LRAT','HAND2','PDGFRB','RBP1','ECM1')

he1=c('AFP','CEBPA', 'CK18','CEBPA', 'CK18', 'HNF4a', 'TAT','ALB', 'CYP3A4','Albumin', 'CALLA', 'CK18', 'CK8','ARG1','Keratin-18', 'Keratin-8','CPS1','CYP2C9','MRP2', 'ZO-1')

cluster = unique(data$cluster)

for (i in cluster){
  gene_set=data[data$cluster==i & data$avg_log2FC>0,'gene']
  gene_intersect=intersect(gene_set,he1)
  if (length(gene_intersect)>0)
  { print(c(i,':',gene_intersect))}
}

# newtype

data=read.csv('New celltype/all_markers_newtype.csv')
p1=jjVolcano(data,  aesCol = c('purple','orange'), size  = 1.5,topGeneN = 10,legend.position = c(0.92,0.93))



# go anlayze

library("clusterProfiler")
library(org.Hs.eg.db)

cluster1=c('Hepatocytes_1','Hepatocytes_2')
cluster2 = cluster[3:7]
for (i in cluster1)
{
  gene_set=data[data$cluster==i & data$avg_log2FC>0,'gene']
  gene_id=bitr(gene_set,fromType='SYMBOL',toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  print(paste0(i,':',length(gene_id$ENTREZID)))
  eG <- enrichGO(gene_id$ENTREZID,OrgDb = org.Hs.eg.db, pvalueCutoff =0.01,qvalueCutoff = 0.01,ont="all",readable =T)
  write.csv(eG,file = paste0(i,'_eG.csv'),row.names = F)
}


library(ggplot2)
library(tidyverse)
filename
h1_eG = read.csv('Hepatocytes_1_eG.csv')


h1_eG <- separate(data=h1_eG, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") 
h1_eG <- separate(data=h1_eG, col=BgRatio, into = c("BR1", "BR2"), sep = "/")
h1_eG <- mutate(h1_eG, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2)))

eGo=h1_eG


eGoBP <- eGo %>%
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>%
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>%
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)



p <- ggplot(eGo10 ,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
p + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')



# find key genes from go analyze and different genes of Hep1 between  HTN and normal
fatty_process=h1_eG[h1_eG$ID=='GO:0006631','geneID']
fatty_process_gene = strsplit(fatty_process,'/')

small_molecule_process=h1_eG[h1_eG$ID=='GO:0044282','geneID']

small_molecule_gene =strsplit(small_molecule_process,'/')[[1]]


differ_gene_Hep1= read.csv('Hepatocytes_1_grup_differential_gene.csv')
up_gene = differ_gene_Hep1[differ_gene_Hep1$p_val<0.005&differ_gene_Hep1$avg_log2FC>0,'X']


fatty_key_gene = intersect(fatty_process_gene[[1]],up_gene)
small_key_gene = intersect(small_molecule_gene,up_gene)


Hep1_go_Data = read.csv('Hepatocytes_1_eG.csv')
new_data  = Hep1_go_Data[,c(2,1,6)]
new_data$all =0
new_data$count =Hep1_go_Data$Count

new_data[new_data$ONTOLOGY=="BP",'all'] = sum(new_data[new_data$ONTOLOGY=="BP",'count'])
new_data[new_data$ONTOLOGY=="MF",'all'] = sum(new_data[new_data$ONTOLOGY=="MF",'count'])

write.csv(new_data,'go.tsv',sep = '\t',row.names = F)

#plot VN picture
a <-list('up genes'= up_gene,
         'fatty gene'= fatty_process_gene[[1]],
         'small_molecule_gene'= small_molecule_gene)
library(gplots)
library(VennDiagram)
library(grid)
library(futile.logger)
library('RColorBrewer')
color <- brewer.pal(3, "Set3")

venn.diagram(
  x = a,
  category.names = c("UP genes" , "fatty genes " , "small molecule genes"),
  filename = 'venn.tiff',
  output=TRUE,  
  imagetype="tiff" , 
  #height = 1000 ,   
  #width = 1000 , 
  resolution = 600, 
  compression = "lzw",   
  lwd = 2, 
  lty = 2, 
  fill = color, 
  col = c("red", 'green', 'blue'), 
  cex = 1.5,  
  fontfamily = "serif",   
  cat.fontface = "bold",
  cat.default.pos = "outer",  
  cat.pos = c(-5, 5, 180), 
  cat.dist = c(0.05, 0.05, 0.05), 
  cat.fontfamily = "serif",
  rotation = 1
)
