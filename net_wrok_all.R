HTN_pro = read.table('HTN_Hep1_promoter.bed',sep ='\t')
HTN_gene = unique(HTN_pro$V5)

nor_pro = read.table('normal_Hep1_promoter.bed',sep ='\t')
nor_gene = unique(nor_pro$V5)



HTN_enh = read.table('HTN_Hep1_enhancer.bed' ,sep = '\t')
nor_enh = read.table('normal_Hep1_enhancer.bed', sep = '\t')


peak_ID = unique(c(HTN_enh$V4, HTN_pro$V4,nor_enh$V4, nor_pro$V4))

motif = read.csv('motif.csv')

motif = motif[motif$X %in% peak_ID, ]


normal_enhancer=read.csv('normal_enhancer_motif.csv')

All_tfs = colnames(motif)

All_tfs=sapply(strsplit(All_tfs,'_'), function(x) x[1])
All_tfs=substring(All_tfs[-1],9)

motif= motif[,-1]
colnames(motif)<- c('ID',All_tfs)


hep_diff_gene = read.csv('Hepatocytes_1_grup_differential_gene.csv')
up_genes = hep_diff_gene[hep_diff_gene$avg_log2FC>0,'X']
down_gens = hep_diff_gene[hep_diff_gene$avg_log2FC<0,'X']


all_RNA_genes = read.csv('integrate_RNA_AverageExpression.csv')

used_TFs = intersect(all_RNA_genes$X,All_tfs)

used_motif = motif[,c('ID',used_TFs)]



nor_TF_enh=NULL

for(ID in nor_enh$V4){
  TF_ID=which(used_motif[used_motif$ID==ID,]==TRUE)
  if (length(TF_ID)>0)
  {
    #TF_names= sapply(TF_ID, function(x) used_TFs[x-1])
    TF_names=sapply(TF_ID, function(x) cbind(used_TFs[x-1],paste0('E_',ID)))
    nor_TF_enh = rbind( nor_TF_enh,t(TF_names))
  }
}


HTN_TF_enh=NULL

for(ID in HTN_enh$V4){
  TF_ID=which(used_motif[used_motif$ID==ID,]==TRUE)
  if (length(TF_ID)>0)
  {
    #TF_names= sapply(TF_ID, function(x) used_TFs[x-1])
    TF_names=sapply(TF_ID, function(x) cbind(used_TFs[x-1],paste0('E_',ID)))
    HTN_TF_enh = rbind( HTN_TF_enh,t(TF_names))
  }
}



used_genes = c(intersect(HTN_gene, nor_gene),intersect(hep_diff_gene$X, c(HTN_gene,nor_gene)))



used_HTN_pro_ID =NULL
used_nor_pro_ID =NULL

for (gene in used_genes)
{
  used_HTN_pro_ID=rbind(used_HTN_pro_ID, HTN_pro[HTN_pro$V5==gene,c('V4','V5')])
  used_nor_pro_ID=rbind(used_nor_pro_ID, nor_pro[nor_pro$V5==gene,c('V4','V5')])
  
}

used_nor_pro_ID =used_nor_pro_ID[!is.na(used_nor_pro_ID$V4),]

used_pro_ID = rbind(used_HTN_pro_ID,used_nor_pro_ID)
colnames(used_pro_ID)<-c('ID','Gene')

use_pro_TF = NULL

for(ID in used_pro_ID$ID){
  TF_ID=which(used_motif[used_motif$ID==ID,]==TRUE)
  if (length(TF_ID)>0)
  {
    #TF_names= sapply(TF_ID, function(x) used_TFs[x-1])
    TF_names=sapply(TF_ID, function(x) cbind(used_TFs[x-1],ID))
    use_pro_TF = rbind( use_pro_TF,t(TF_names))
  }
}



library(dplyr)
colnames(use_pro_TF)<-c('TF','ID')
use_pro_TF = as.data.frame(use_pro_TF)
use_pro_TF$ID = as.integer(use_pro_TF$ID) 

pro_ID_TF =  inner_join(used_pro_ID,use_pro_TF,by="ID")


HTN_TF_enh =as.data.frame(HTN_TF_enh)
HTN_TF_enh$enhancer='HTN_E'

nor_TF_enh =as.data.frame(nor_TF_enh)
nor_TF_enh$enhancer='Nor_E'

TF_en = rbind(nor_TF_enh,HTN_TF_enh)

colnames(TF_en)<-c('TF','E_ID','enhancer')

#pro_ID_TF_enh =  inner_join(pro_ID_TF,TF_en,by="TF")

pro_ID_TF_enh =  inner_join(pro_ID_TF,TF_en[,c('TF','enhancer')],by="TF")

duplicates <- duplicated(pro_ID_TF_enh)

pro_ID_TF_enh_no_duplicates <- pro_ID_TF_enh[!duplicates, ]


sample_gene =NULL

HTN = rbind('HTN',used_genes[used_genes%in%HTN_gene])
Nor = rbind('Normal',used_genes[used_genes%in%nor_gene])

sample_gene = as.data.frame(rbind(t(HTN),t(Nor))) 
colnames(sample_gene)<-c('sample','Gene')

network = inner_join(sample_gene,pro_ID_TF_enh_no_duplicates,by="Gene")

network_data = network[ ,c("sample","Gene","TF","enhancer")]

network_df <- to_lodes_form(network_data[,1:ncol(network_data)],
                    axes = 1:ncol(network_data),
                    id = "value")


library(ggplot2)
library(ggalluvial)

write.table(network_data,'network_data.csv',row.names = F,sep = '\t')

col<- rep(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5"), 50)

pdf("network.pdf",width = 10, height = 20)

ggplot(network_df, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+
  geom_flow(width = 0.3,
            curve_type = "sine",
            alpha = 0.5,
            color = 'white',
            size = 0.1)+
  geom_stratum(width = 0.28)+
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col )+
  theme_void()+
  theme(legend.position = 'none')

dev.off()


# plot single network

library(ggraph)
library(tidygraph)
library(igraph)
library(ggplot2)

CYP3A4_net = NULL


CYP3A4_TF_ID = which(used_motif[used_motif$ID==70606,]==TRUE)-1

CYP3A4_TF=  sapply(CYP3A4_TF_ID, function(x) cbind(used_TFs[x],'CYP3A4'))

CYP3A4_net= as.data.frame(t(CYP3A4_TF))
colnames(CYP3A4_net)<-c('source','target')

CYP3A4_inf = as.data.frame(cbind('CYP3A4','gene'))

CYP3A4_inf= rbind(CYP3A4_inf,t(rbind(CYP3A4_net$source,'TF')))

for (TF in CYP3A4_TF_ID)
{
  TF_peaK_ID = which(used_motif[,TF+1]==TRUE)
  en_ID = intersect(TF_peaK_ID,HTN_enh$V4)
  if (length(en_ID)>0)
  {
    net=as.data.frame(cbind(en_ID,used_TFs[TF])) 
    colnames(net)<-c('source','target')
    inf = t(rbind(net$source,'enhancer'))
    CYP3A4_net= rbind(CYP3A4_net,net)
    CYP3A4_inf=rbind(CYP3A4_inf,inf)
  }
 
}

CYP3A4_inf= as.data.frame(CYP3A4_inf)

colnames(CYP3A4_inf)<-c('name','class')

g <- tbl_graph(nodes = CYP3A4_inf,edges = CYP3A4_net,directed = T)

V(g)$degree <- degree(g)

mycol <- c("#FF8901","#00C5FF","#FF5485")



p6 <- ggraph(g, layout = 'linear',circular = TRUE)+
  geom_edge_arc(colour="grey50",width=1,alpha=0.3)+
  geom_node_point(aes(color=class),size=6,alpha=0.8)+
  scale_colour_manual(values = mycol)+
  geom_node_text(aes(x = 1.05 * x,y = 1.05 * y,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
                     label=name),
                 nudge_y = 0,
                 hjust = 'outward',
                 repel = F,
                 size=2.5)+
  coord_fixed(clip = "off")+
  theme_graph()


WT1_peak_ID  = which(used_motif[,'WT1']==TRUE)

EGR3_peak_ID = used_motif[which(used_motif[,'EGR3']==TRUE),'ID']

EGR3_nor_gene= paste0('N_',nor_pro[nor_pro$V4 %in% EGR3_peak_ID,'V5'])
EGR3_HTN_gene =paste0('H_',HTN_pro[HTN_pro$V4 %in% EGR3_peak_ID,'V5']) 

EGR3_nor_enh =nor_enh[nor_enh$V4 %in% EGR3_peak_ID,'V4']
EGR3_HTN_enh =HTN_enh[HTN_enh$V4 %in% EGR3_peak_ID,'V4']

EGR3_net = as.data.frame(t(rbind('EGR3',c(EGR3_nor_gene,EGR3_HTN_gene))))

colnames(EGR3_net)<-c('source','target')

net = as.data.frame(t(rbind(EGR3_HTN_enh,'EGR3')))

colnames(net)<-c('source','target')
EGR3_net =rbind(EGR3_net,net)

EGR3_duplicates <- duplicated(EGR3_net)
EGR3_net1 =EGR3_net[!EGR3_duplicates,]




EGR3_inf = as.data.frame(cbind('EGR3','TF'))
colnames(EGR3_inf)<-c('name','class')
inf = as.data.frame(t(rbind(net$source,'HTN_enhancer')))
colnames(inf)<-c('name','class')
EGR3_inf = rbind(EGR3_inf,inf)

inf = as.data.frame(t(rbind(c(EGR3_nor_gene,EGR3_HTN_gene),'gene')))
colnames(inf)<-c('name','class')
EGR3_inf = rbind(EGR3_inf,inf)
EGR3_inf1 =EGR3_inf[!duplicated(EGR3_inf),]



g <- tbl_graph(nodes = EGR3_inf1,edges = EGR3_net1,directed = F)
V(g)$degree <- degree(g)


ggraph(g, layout = 'linear',circular = TRUE)+
  geom_edge_arc(colour="grey50",width=1,alpha=0.3)+
  geom_node_point(aes(color=class),size=6,alpha=0.8)+
  scale_colour_manual(values = mycol)+
  geom_node_text(aes(x = 1.05 * x,y = 1.05 * y,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
                     label=name),
                 nudge_y = 0,
                 hjust = 'outward',
                 repel = F,
                 size=2.5)+
  coord_fixed(clip = "off")+
  theme_graph()
