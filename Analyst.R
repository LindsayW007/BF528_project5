library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

pbmc <- readRDS('pbmc_programmer.rda')

#Identify markers for each 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers_sig <- pbmc.markers[pbmc.markers$p_val_adj<0.05, ]
top10_marker <- pbmc.markers_sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(pbmc.markers_sig, file="gene_marker_sig.csv")

###Marker Genes for Cluster Labeling
#Alpha: GCG
#Beta: INS, MAFA
#Gamma: PPY
#Delta: SST
#Endothelial: SPARC, VWF
#Acinar: PRSS1, CPA1
#Ductal: KRT19, SOX4(ref)
#Quiescent Stellate: RGS5
#Activated Stellate: PDGFRA
#Schwann cell: SOX10
#Macrophage: SDS

###
#--------------------------------------------------------------------------
#Cluster  |           0|           1|    2|     3|     4|      5|      6|
#--------------------------------------------------------------------------
#Cell Type| Delta/Gamma| Alpha/Gamma| Beta| Gamma| Alpha| Acinar| Ductal| 
#--------------------------------------------------------------------------
#Cluster  |          7|    8|      9|          10|          11|         12|
#--------------------------------------------------------------------------
#Cell Type|Endothelial| Beta| Ductal| Delta/Gamma| Endothelial| Macrophage|
#--------------------------------------
###

gene_markers <- c('GCG', 'INS', 'MAFA', 'PPY', 'SST',
                  'SPARC', 'VWF', 'PRSS1', 'CPA1', 'SOX4', 'SDS')
cluster_label <- c('D/G_1', 'A/G', 'Beta_1', 'Gamma', 'Alpha', 'Acinar', 'Ductal_1',
                   'Endothelial_1', 'Beta_2', 'Ductal_2', 'D/G_2', 'Endothelial_2', 'Macrophage')

sig_genes <- c("REG1B", "PLVAP", "PRSS3P2", "CTRB2", "COL1A1", "ALB", "PRSS1", "REG1A", "CTRB1", "IGFBP5")

names(cluster_label) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cluster_label)
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
png(filename = 'umap_plot.png', width=800, height=600)
umap_plot
dev.off()

#Heatmap for marker genes from literature
heatmap_lit_marker<- DoHeatmap(pbmc, features = gene_markers, size=0, draw.lines=FALSE)
png(filename = 'heatmap_lit_marker.png', width=1000, height=500)
heatmap_lit_marker
dev.off()

#Violn Plot for marker genes from literature 
vln_plot_lit_marker <- VlnPlot(pbmc, features = gene_markers, pt.size=0, combine=TRUE)
png(filename = 'vln_plot_lit_marker.png', width=2000, height=800)
vln_plot_lit_marker
dev.off()


###############Identify Novel Marker Genes
marker_alpha <- subset(top10_marker$gene, top10_marker$cluster==4)
marker_alpha <- marker_alpha[marker_alpha != 'GCG']
marker_alpha <- c('GCG', marker_alpha)
nm_alpha <- DoHeatmap(pbmc, features = marker_alpha, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Alpha Cells')
png(filename = 'nm_alpha.png', width=1000, height=500)
nm_alpha
dev.off()


marker_beta <- subset(top10_marker$gene, top10_marker$cluster==2 | top10_marker$cluster==8)
marker_beta <- marker_beta[marker_beta != 'INS' & marker_beta != 'MAFA']
marker_beta <- c('INS', 'MAFA', marker_beta)
nm_beta <- DoHeatmap(pbmc, features = marker_beta, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Beta Cells')
png(filename = 'nm_beta.png', width=1000, height=500)
nm_beta
dev.off()


marker_gamma <- subset(top10_marker$gene, top10_marker$cluster==3)
marker_gamma <- marker_gamma[marker_gamma != 'PPY']
marker_gamma <- c('PPY', marker_gamma)
nm_gamma <-DoHeatmap(pbmc, features = marker_gamma, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Gamma Cells')
png(filename = 'nm_gamma.png', width=1000, height=500)
nm_gamma
dev.off()


marker_delta <- subset(top10_marker$gene, top10_marker$cluster==0 | top10_marker$cluster==10)
marker_delta <- marker_delta[marker_delta != 'SST']
marker_delta <- c('SST', marker_delta)
nm_delta <-DoHeatmap(pbmc, features = marker_delta, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Delta Cells')
png(filename = 'nm_delta.png', width=1000, height=500)
nm_delta
dev.off()


marker_acinar <- subset(top10_marker$gene, top10_marker$cluster==5)
marker_acinar <- marker_acinar[marker_acinar != 'PRSS1' & marker_acinar != 'CPA1']
marker_acinar <- c('PRSS1', 'CPA1', marker_acinar)
nm_acinar <-DoHeatmap(pbmc, features = marker_acinar, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Acinar Cells')
png(filename = 'nm_acinar.png', width=1000, height=500)
nm_acinar
dev.off()


marker_ductal <- subset(top10_marker$gene, top10_marker$cluster==6 | top10_marker$cluster==9)
marker_ductal <- marker_ductal[marker_ductal != 'KRT19' & marker_ductal != 'SOX4']
marker_ductal <- c('KRT19', 'SOX4', marker_ductal)
nm_ductal<-DoHeatmap(pbmc, features = marker_ductal, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Ductal Cells')
png(filename = 'nm_ductal.png', width=1000, height=500)
nm_ductal
dev.off()


marker_endothelial <- subset(top10_marker$gene, top10_marker$cluster==7 | top10_marker$cluster==11)
marker_endothelial <- marker_endothelial[marker_endothelial != 'SPARC' & marker_endothelial != 'VWF']
marker_endothelial <- c('SPARC', 'VWF', marker_endothelial)
nm_endothelial <-DoHeatmap(pbmc, features = marker_endothelial, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Endothelial Cells')
png(filename = 'nm_endothelial.png', width=1000, height=500)
nm_endothelial
dev.off()


marker_macrophage <- subset(top10_marker$gene, top10_marker$cluster==12)
marker_macrophage <- marker_macrophage[marker_macrophage != 'SDS']
marker_macrophage <- c('SDS', marker_macrophage)
nm_macrophage <-DoHeatmap(pbmc, features = marker_macrophage, size=0, draw.lines=FALSE) + labs(title='Novel Marker for Macrophage Cells')
png(filename = 'nm_macrophagel.png', width=1000, height=500)
nm_macrophage
dev.off()

#Create table for potential novel markers
m_alpha <- subset(top10_marker, top10_marker$cluster==4)
m_alpha <- subset(m_alpha, m_alpha$gene!= 'GCG')
m_alpha <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==4 & pbmc.markers_sig$gene=='GCG'), m_alpha)
m_alpha$marker_for <- 'alpha'

m_beta <- subset(top10_marker, top10_marker$cluster==2 | top10_marker$cluster==8)
m_beta <- subset(m_beta, m_beta$gene!= 'INS' & m_beta$gene!= 'MAFA')
m_beta <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==2 & pbmc.markers_sig$gene=='INS'), 
                subset(pbmc.markers_sig, pbmc.markers_sig$cluster==8 & pbmc.markers_sig$gene=='MAFA'), m_beta)
m_beta$marker_for <- 'beta'

m_gamma <- subset(top10_marker, top10_marker$cluster==3)
m_gamma <- subset(m_gamma, m_gamma$gene!= 'PPY')
m_gamma <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==3 & pbmc.markers_sig$gene=='PPY'), m_gamma)
m_gamma$marker_for <- 'gamma'

m_delta <- subset(top10_marker, top10_marker$cluster==0 | top10_marker$cluster==10)
m_delta <- subset(m_delta, m_delta$gene!= 'SST')
m_delta <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==0 & pbmc.markers_sig$gene=='SST'), 
                 subset(pbmc.markers_sig, pbmc.markers_sig$cluster==10 & pbmc.markers_sig$gene=='SST'), m_delta)
m_delta$marker_for <- 'delta'

m_acinar <- subset(top10_marker, top10_marker$cluster==5)
m_acinar <- subset(m_acinar, m_acinar$gene!= 'PRSS1' & m_acinar$gene!= 'CPA1')
m_acinar <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==5 & pbmc.markers_sig$gene=='PRSS1'), 
                subset(pbmc.markers_sig, pbmc.markers_sig$cluster==5 & pbmc.markers_sig$gene=='CPA1'), m_acinar)
m_acinar$marker_for <- 'acinar'

m_endothelial <- subset(top10_marker, top10_marker$cluster==7 | top10_marker$cluster==11)
m_endothelial <- subset(m_endothelial, m_endothelial$gene!= 'SPARC' & m_endothelial$gene!= 'VWF')
m_endothelial <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==7 & pbmc.markers_sig$gene=='SPARC'), 
                 subset(pbmc.markers_sig, pbmc.markers_sig$cluster==11 & pbmc.markers_sig$gene=='SPARC'), 
                 subset(pbmc.markers_sig, pbmc.markers_sig$cluster==11 & pbmc.markers_sig$gene=='VWF'), 
                 m_endothelial)
m_endothelial$marker_for <- 'endothelial'

m_macrophage <- subset(top10_marker, top10_marker$cluster==12)
m_macrophage <- subset(m_macrophage, m_macrophage$gene!= 'SDS')
m_macrophage <- rbind(subset(pbmc.markers_sig, pbmc.markers_sig$cluster==12 & pbmc.markers_sig$gene=='SDS'), m_macrophage)
m_macrophage$marker_for <- 'macrophage'


#Write the file with potential novel markers 
novel_markers <- rbind(m_alpha, m_beta, m_delta, m_gamma, m_acinar, m_endothelial, m_macrophage)
novel_markers <- subset(novel_markers, select = -c(cluster))
write.csv(novel_markers, 'novel_markers.csv', row.names=FALSE)
          