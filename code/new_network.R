diff_gene = read.csv('./RNA/Hepatocytes_1_grup_differential_gene.csv')
HTN_pro = read.table('./promoter_enhancer/HTN_Hep1_promoter.bed',sep = '\t')
HTN_enh = read.table('./promoter_enhancer/HTN_Hep1_enhancer.bed',sep = '\t')
HTN_GeneWithPromoter = unique(HTN_pro$V5)
up_gene = diff_gene[diff_gene$avg_log2FC>0,'X']
used_gene = intersect(up_gene,HTN_GeneWithPromoter)

gene_promterID= HTN_pro[HTN_pro$V5%in%used_gene,c('V4','V5')]
colnames(gene_promterID)<-c('ID','Gene')

motif = read.csv('all_motif.csv')
motif = motif[motif$X %in% c(HTN_enh$V4, HTN_pro$V4),]
All_tfs = colnames(motif)
All_tfs=sapply(strsplit(All_tfs,'_'), function(x) x[1])
All_tfs=substring(All_tfs[-1],9)
colnames(motif)<- c('ID',All_tfs)

use_pro_TF =''
for(ID in gene_promterID$ID){
  TF_ID=which(motif[motif$ID==ID,]==TRUE)
  if (length(TF_ID)>0)
  {
    TF_names=sapply(TF_ID, function(x) cbind(All_tfs[x-1],ID))
    use_pro_TF = rbind( use_pro_TF,t(TF_names))
  }
}

colnames(use_pro_TF)<-c('TF','ID')
use_pro_TF = as.data.frame(use_pro_TF[-1,])
use_pro_TF$ID = as.integer(use_pro_TF$ID) 
library(dplyr)
pro_ID_TF =  inner_join(gene_promterID,use_pro_TF,by="ID")

all_RNA_genes = read.csv('./RNA/integrate_RNA_AverageExpression.csv')

pro_ID_TF = pro_ID_TF[pro_ID_TF$TF%in%all_RNA_genes$X,]
write.csv(pro_ID_TF,'promoter_net.csv',row.names = F)


#__________enhancer_Tfs______


library(readxl)
library(GenomicRanges)
data <- read_excel("GeneHancer_Version_4-4.xlsx", sheet = "Sheet1")

HTN_enhancer_Grange=GRanges(
  seqnames = HTN_enh$V1,
  ranges = IRanges(start = HTN_enh$V2, end = HTN_enh$V3)
)

all_gene_card_enhancer_Grange=GRanges(
  seqnames = data$chrom,
  ranges = IRanges(start = data$start, end = data$end)
)

all_overlop = findOverlaps(all_gene_card_enhancer_Grange,HTN_enhancer_Grange)

geneWithEnhancer = NULL
for (gene in up_gene)
{
  genecard_enhancer=data[which(grepl(gene, data$attributes)),] 
  genecard_enhancer_Grange = GRanges(
    seqnames = genecard_enhancer$chrom,
    ranges = IRanges(start = genecard_enhancer$start, end = genecard_enhancer$end)
  )
  overlop = findOverlaps(genecard_enhancer_Grange,HTN_enhancer_Grange)
  if (length(overlop)>1)
  {
    gene_card_ID =sapply(genecard_enhancer[overlop@from,'attributes'], function(x) substring(x, 15, 25))
    ID = HTN_enh[overlop@to,'V4']
    geneWithEnhancer=rbind(geneWithEnhancer,cbind(gene,ID,gene_card_ID))
    #break
  }
}

geneWithEnhancer =as.data.frame(geneWithEnhancer)
colnames(geneWithEnhancer)<-c("gene","ID","genehancer_id")
geneWithEnhancer$ID = as.integer(geneWithEnhancer$ID)


enhancer_tf = NULL

for(ID in geneWithEnhancer$ID){
  TF_ID=which(motif[motif$ID==ID,]==TRUE)
  if (length(TF_ID)>0)
  {
    TF_names=sapply(TF_ID, function(x) cbind(All_tfs[x-1],ID))
    enhancer_tf = rbind(enhancer_tf,t(TF_names))
  }
}
enhancer_tf =as.data.frame(enhancer_tf)
colnames(enhancer_tf)<-c("TF","ID")
enhancer_tf$ID = as.integer(enhancer_tf$ID)
gene_enhancer_ID_TF = inner_join(geneWithEnhancer,enhancer_tf,by="ID")
gene_enhancer_ID_TF = gene_enhancer_ID_TF[gene_enhancer_ID_TF$TF%in%all_RNA_genes$X,]


Location=HTN_enh[HTN_enh$V4%in%gene_enhancer_ID_TF$ID,c('V1','V2','V3','V4')] 
write.csv(Location,'enhancer_location.csv',row.names = F)

write.csv(gene_enhancer_ID_TF,'enhancer_net.csv',row.names = F)


#--------------plot network for promoter---------------

library(ggplot2)
library(ggalluvial)

promoter_net = pro_ID_TF[ ,c("Gene","ID","TF")]
promoter_location = HTN_pro[HTN_pro$V4 %in% promoter_net$ID,c('V1','V2','V3','V4','V5')]


  
promoter_net_df <- to_lodes_form(promoter_net[,1:ncol(promoter_net)],
                            axes = 1:ncol(promoter_net),
                            id = "value")
col<- rep(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5"), 50)

pdf("promoter.pdf",width = 5, height = 5)




ggplot(promoter_net_df, aes(x = x, fill = stratum, label = stratum, 
                            stratum = stratum, alluvium = value)) +
  geom_flow(width = 0.3,  # Adjusted width to make columns wider
            curve_type = "sine", 
            alpha = 0.5, 
            color = 'white') +
  geom_stratum(width = 0.28) +  # Adjusted width to make columns wider
  geom_text(stat = 'stratum', size = 4, color = 'black') +  # Increase text size
  scale_fill_manual(values = col) +
  theme_void() +
  theme(legend.position = 'none')

dev.off()

#-------plot enhancer net---------

gene_enhancer_ID_TF$id=paste0(gene_enhancer_ID_TF$ID,'(',gene_enhancer_ID_TF$genehancer_id,')')
gene_net = gene_enhancer_ID_TF[ ,c("gene","ID","TF")]



enhancer_net_df <- to_lodes_form(gene_net[,1:ncol(gene_net)],
                                 axes = 1:ncol(gene_net),
                                 id = "value")
col<- rep(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5"), 50)
pdf("enhancer1.pdf",width = 5, height = 10)
ggplot(enhancer_net_df, aes(x = x, fill = stratum, label = stratum, 
                            stratum = stratum, alluvium = value)) +
  geom_flow(width = 0.4,  # Adjusted width to make columns wider
            curve_type = "sine", 
            alpha = 0.5, 
            color = 'white') +
  geom_stratum(width = 0.45) +  # Adjusted width to make columns wider
  geom_text(stat = 'stratum', size = 3.5, color = 'black') +  # Increase text size
  scale_fill_manual(values = col) +
  theme_void() +
  theme(legend.position = 'none')

dev.off()



#！！！！！！！！！！chip-seq--------------



liver_chip = read.csv('liver_chip_seq.tsv',header = T, sep = '\t')
HepG2_chip = read.csv('HepG2_chip_seq.tsv', header = T, sep='\t')

all_chip = rbind(liver_chip,HepG2_chip)
all_chip_1 = all_chip[c("File.download.URL","Experiment.target","Biosample.term.name","File.format","File.accession")]
all_chip_1$Experiment.target=sapply(strsplit(all_chip_1$Experiment.target,'-'), function(x) x[1])
all_chip_1= all_chip_1[all_chip_1$File.format=='bed narrowPeak',]

library(rtracklayer)

Chip_seq= read.csv('./Chip_seq/Oth.Liv.05.FOXP1.AllCell.bed',sep = '\t',header = F,skip = 1)
Chip_seq_Grange=  GRanges(
  seqnames = Chip_seq$V1,
  ranges = IRanges(start = Chip_seq$V2, end = Chip_seq$V3)
)

TF = 'FOXP1'

TF_in_promoter=pro_ID_TF[pro_ID_TF$TF==TF,'ID']

TF_in_promoter_Grange =GRanges(
  seqnames = HTN_pro[HTN_pro$V4%in%TF_in_promoter,'V1'],
  ranges = IRanges(HTN_pro[HTN_pro$V4%in%TF_in_promoter,'V2'], 
                   end = HTN_pro[HTN_pro$V4%in%TF_in_promoter,'V3'])
)

findOverlaps(TF_in_promoter_Grange,Chip_seq_Grange)


