library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library(fishpond)
library(EnsDb.Hsapiens.v79)
library(ggplot2)
library(ggrepel)
library(scales)

#Read alevin output
files = file.path("/projectnb/bf528/users/van-gogh/project_4/data/salmon_output_merged_drop_100/alevin/quants_mat.gz")
txi <- tximport(files, type="alevin")
table <- txi$counts

#Change matrix index from ensembleid to gene name
ensembl_id <-  rownames(table) 
ensembl_id <- sub("[.][0-9]*","",ensembl_id)
rownames(table) <- ensembl_id
#Mapping ensembleid to gene name 
Gene_symbol <- select(EnsDb.Hsapiens.v79, keys=ensembl_id, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
table <- table[rownames(table) %in% Gene_symbol$GENEID,]
rownames(table) <- Gene_symbol$SYMBOL

#Create Seurat Object
pbmc <- CreateSeuratObject(counts = table, min.cells = 3, min.features = 200, project = "project5")
pbmc

#QC for identifying the percentage of reads that map to the mitochondrial genome
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

####Filter out low-quality cells
# Visualize QC metrics as a violin plot
violin_plot_QC <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 'violin_plot_QC.png', width=500, height=600)
violin_plot_QC
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
feature_plot_QC <- plot1 + plot2
png(filename = 'feature_plot_QC.png', width=1000, height=600)
feature_plot_QC
dev.off()

#QC: filter cells that have unique feature counts over 4,000 or less than 200
#QC: filter cells that have >40% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 40)

#Check the distribution after filtering
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

###Filter out low-variance gene
#Normalization 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#Select highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
#Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
sig_genes_QC <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(filename = 'sig_genes_QC.png', width=1000, height=600)
sig_genes_QC
dev.off()

#Scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA for scaled data
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")

#Determine dimentionality of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
elbow_QC <- ElbowPlot(pbmc)
png(filename = 'elbow_QC.png', width=600, height=600)
elbow_QC
dev.off()

#cluster the cells with 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

#UMAP visulaization
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
#saveRDS(pbmc, file="pbmc_programmer.rda")

#Cell count per cluster summary
clusters <- data.frame(table(Idents(pbmc)))
names(clusters)[names(clusters) == 'Var1'] <- 'cluster'
clusters$percentage <- clusters$Freq/sum(clusters$Freq)
clusters$percentage <- label_percent()(clusters$percentage)
clusters['percentage'] <- round(clusters['percentage'], 2)
clusters$pos <- (sum(clusters$Freq)-(cumsum(c(0, clusters$Freq)) + c(clusters$Freq / 2, .01)))[1:nrow(clusters)]

bar <-  ggplot(clusters, aes(x = cluster, y = Freq, fill = cluster)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = 'Cell Count per Cluster') +
  theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey", size=0.1), plot.title = element_text(size=18))
png(filename = 'bar_cluster.png', width=700, height=700)
bar
dev.off()

pie <- ggplot(clusters, aes(x = '', y = Freq, fill = cluster)) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
  geom_text_repel(aes(x = 1.4, y = pos, label = percentage), 
                  nudge_x = .3, 
                  segment.size = .3, 
                  show.legend = FALSE) +
  labs(title = 'Cell Count Percentage per Cluster') + 
  theme_void() +
  theme(plot.title = element_text(size=18))
png(filename = 'pie_cluster.png', width=700, height=700)
pie
dev.off()          
