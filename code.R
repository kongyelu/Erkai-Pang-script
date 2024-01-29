setwd('')
library(Seurat)
library(tidyverse)
fs=list.files('./','^GSM')
library(stringr)
samples=str_split(fs,'_',simplify = T)[,1]

lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste(str_split(y[1],'_',simplify = T)[,1],
               collapse = '')
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"genes.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})

folders=list.files('./',pattern='[0-9]$')
metad <- read.csv('meta.csv',header = F)
colnames(metad) <- c('GSM','group1','group2')
rownames(metad) <- metad$GSM

for (i in folders) {
  assign(i,CreateSeuratObject(counts = Read10X(i), project = i ,min.cells = 3, min.features = 200))
  assign(i,AddMetaData(object = get(i),metadata = 'GSE220243',col.name = 'Study'))
  assign(i,AddMetaData(object = get(i),metadata = metad[i,'group1'],col.name = 'group1'))
  assign(i,AddMetaData(object = get(i),metadata = metad[i,'group2'],col.name = 'group2'))
}
rm(samples,fs,folders,i,metad)
file <- ls()
gc()

GSE_raw_list <- lapply(file,FUN = get)
set.resolutions <- seq(0.2, 2, by = 0.1)

doubletDetect <- function(Seurat.object, PCs, doublet.rate = 0.076, annotation = "seurat_clusters", pN_value = 0.25, GT = FALSE, sct = FALSE){
  library(DoubletFinder) 
  T1.sample <- paramSweep_v3(Seurat.object, PCs = PCs, sct = sct)
  T1.gt.calls <- summarizeSweep(T1.sample, GT = GT)
  
  
  sweep.T1 <- find.pK(T1.gt.calls)
  pK_value <- as.numeric(as.character(sweep.T1$pK[sweep.T1$BCmetric == max(sweep.T1$BCmetric)])) 
  
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  plot(x = as.numeric(as.character(sweep.T1$pK)), y = sweep.T1$BCmetric, pch = 16,type="b", col = "blue",lty=1)
  abline(v=pK_value,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_value,max(sweep.T1$BCmetric),as.character(pK_value),pos = 4,col = "red")

  nExp_poi <- round(doublet.rate*nrow(Seurat.object@meta.data))
  if(is.null(annotation)){
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
    label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi)
  }else{
    annotations <- Seurat.object@meta.data[, annotation]
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    pANN_value <- paste0("pANN_", pN_value, "_", pK_value, '_', nExp_poi)
    
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = sct)
    label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi.adj)
  }
  Seurat.object@meta.data$Doublet <- Seurat.object@meta.data[, label]
  
  return(Seurat.object)
}
set.resolutions <- seq(0.2, 2, by = 0.1)
preprocess <- function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- subset(x,  subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
  x <- SCTransform(x, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = TRUE,vst.flavor = "v2")
  x <- RunPCA(x, npcs = 100, verbose = TRUE)
  x  <- FindNeighbors(object = x , dims = 1:50, verbose = TRUE)
  x  <- FindClusters(object = x , resolution = set.resolutions, verbose = TRUE) 
  x  <- RunUMAP(x , dims = 1:50)
  x
}
GSE_raw_list <- lapply(GSE_raw_list, preprocess)
GSE_raw_list <- lapply(GSE_raw_list,doubletDetect,PCs = 1:50, annotation = "SCT_snn_res.0.5", sct = T)
data_merge=merge(x=GSE_raw_list[[1]],y=GSE_raw_list[2:length(GSE_raw_list)])
data_merge_remove_doublet <- subset(data_merge, subset = Doublet == "Singlet")

DefaultAssay(data_merge_remove_doublet) <- "RNA"
data_merge_remove_doublet.list <- SplitObject(data_merge_remove_doublet, split.by = "orig.ident")
data_merge_remove_doublet.list=data_merge_remove_doublet.list[which(lapply(data_merge_remove_doublet.list,FUN = function(x){
  dim(x@meta.data)[1]
})>30)]
data_merge_remove_doublet.list <- lapply(X = data_merge_remove_doublet.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T,vst.flavor = "v2")
})
data_merge_remove_doublet.scRNA <- merge(data_merge_remove_doublet.list[[1]], y = data_merge_remove_doublet.list[2:length(data_merge_remove_doublet.list)], project = "data_merge_remove_doublet")
DefaultAssay(data_merge_remove_doublet.scRNA) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = data_merge_remove_doublet.list)
VariableFeatures(data_merge_remove_doublet.scRNA) <- seurat.features.SCT

Harmony.integration.reduceDimension <- function(seurat.object, set.resolutions, assay = "RNA", nfeatures = 3000, PC = 50, npcs = 100){
  require(Seurat)
  require(harmony)
  require(dplyr)
  DefaultAssay(seurat.object) <- assay
  
  seurat.object <- RunPCA(seurat.object, verbose = FALSE, npcs = npcs)
  p <- ElbowPlot(object = seurat.object, ndims = npcs)
  
  seurat.object <- RunHarmony(object = seurat.object, group.by.vars = "orig.ident", assay.use=assay, verbose = FALSE)
  seurat.object <- RunUMAP(seurat.object, reduction = "harmony", dims = 1:PC, verbose = FALSE)
  seurat.object <- RunTSNE(seurat.object, reduction = "harmony", dims = 1:PC, verbose = FALSE)
  seurat.object <- FindNeighbors(seurat.object, dims = 1:PC, reduction = "harmony", verbose = FALSE) 
  seurat.object <- FindClusters(seurat.object, resolution = set.resolutions, verbose = FALSE)
  return(seurat.object)
}

data.harmony <- Harmony.integration.reduceDimension(seurat.object = data_merge_remove_doublet.scRNA, assay = "SCT", 
                                                    set.resolutions = seq(0.2, 2, by = 0.1), PC = 30, npcs = 50)
save(data.harmony,file='harmony.RData')

setwd('')
library(Seurat)
load('harmony.RData')
features <- c('CHI3L1','CHI3L2','CHI3L1','CHI3L2','MT1G','MT1H','MT1E','MT1X','CYTL1','FRZB','CHAD','CLEC3A','S100A1','S100B','PRG4','SEMA3A',
              'OGN','COL2A1','COL3A1','IGFBP5','S100A4','MSMP','COL1A1','COL1A2','IGFBP5','S100A4','PRG4','SEMA3A','SOX9','COL9A3','COL10A1',
              'IBSP','SOX9','JUN','MMP3','RGS16','CILP','CILP2','COL2A1')
features <- intersect(features,rownames(data.harmony))

FeaturePlot(data.harmony, features = 'GNLY', cols = c("lightgrey", "red"), reduction = 'umap',raster=FALSE)
DimPlot(data.harmony,group.by = 'SCT_snn_res.1.5',raster=FALSE)
pdf('m.pdf')
res <- sapply(features, function(x) {
  p <- FeaturePlot(data.harmony, features = x, cols = c("lightgrey", "red"), reduction = 'umap',raster=FALSE)
  print(p)

})
dev.off()
set.resolutions <- seq(0.2, 2, by = 0.1)
pdf('d.pdf')
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.harmony, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x),raster=F)
  print(p)
})
dev.off()

data.harmony@meta.data$seurat_clusters=data.harmony@meta.data$SCT_snn_res.0.8
cluster.label <- data.harmony@meta.data$seurat_clusters
cluster.label <- gsub("^2$", "Regulatory Chondrocytes-1", cluster.label)

cluster.label <- gsub("^1$", "Regulatory Chondrocytes-2", cluster.label)
cluster.label <- gsub("^17$", "Regulatory Chondrocytes-2", cluster.label)
cluster.label <- gsub("^18$", "Regulatory Chondrocytes-2", cluster.label)

cluster.label <- gsub("^0$", "Effector Chondrocytes", cluster.label)
cluster.label <- gsub("^5$", "Effector Chondrocytes", cluster.label)
cluster.label <- gsub("^6$", "Effector Chondrocytes", cluster.label)
cluster.label <- gsub("^13$", "Effector Chondrocytes", cluster.label)
cluster.label <- gsub("^15$", "Effector Chondrocytes", cluster.label)
cluster.label <- gsub("^116$", "Effector Chondrocytes", cluster.label)

cluster.label <- gsub("^4$", "Pre-fibrocartilage Chondrocytes", cluster.label)
cluster.label <- gsub("^3$", "Pre-fibrocartilage Chondrocytes", cluster.label)
cluster.label <- gsub("^8$", "Pre-fibrocartilage Chondrocytes", cluster.label)
cluster.label <- gsub("^14$", "Pre-fibrocartilage Chondrocytes", cluster.label)
cluster.label <- gsub("^16$", "Pre-fibrocartilage Chondrocytes", cluster.label)



cluster.label <- gsub("^9$", "Fibrocartilage Chondrocytes", cluster.label)
cluster.label <- gsub("^12$", "Fibrocartilage Chondrocytes", cluster.label)
cluster.label <- gsub("^10$", "Hypertrophic Chondrocytes", cluster.label)

cluster.label <- gsub("^7$", "Homeostatic Chodrocytes", cluster.label)

cluster.label <- gsub("^11$", "Reparative Chondrocytes", cluster.label)

data.harmony <- AddMetaData(data.harmony, cluster.label, col.name = "cellType_low")
table(data.harmony$cellType_low)
selcet_cluster <- setdiff(levels(data.harmony@meta.data$seurat_clusters),c())
data_harmony_pro <- subset(data.harmony, subset = seurat_clusters %in% selcet_cluster)
table(data_harmony_pro$cellType_low)
data_harmony_pro$seurat_clusters <- factor(data_harmony_pro$seurat_clusters, levels = selcet_cluster)
data_harmony_pro$cellType_low<- factor(data_harmony_pro$cellType_low,levels = c('Regulatory Chondrocytes-1','Regulatory Chondrocytes-2','Effector Chondrocytes',
                                                                                'Pre-fibrocartilage Chondrocytes','Fibrocartilage Chondrocytes','Pre-hypertrophic Chondrocytes',
                                                                                'Hypertrophic Chondrocytes','Homeostatic Chodrocytes','Reparative Chondrocytes'))
pdf('dimplot.pdf',height = 10,width = 12.5)
DimPlot(data_harmony_pro,group.by ='cellType_low',raster=F,label = T)
dev.off()
save(data_harmony_pro, file = "data_harmony_pro.RData")
library(ggplot2)

pdf('dot.pdf',width = 25,height = 10)
DoHeatmap(data_harmony_pro,  features = features, group.by = "cellType_low")
dev.off()
Idents(data_harmony_pro) <- data_harmony_pro$SCT_snn_res.0.8
pdf('vlnPlot.pdf',width = 20,height = 10)
VlnPlot(data_harmony_pro,features = c("nCount_RNA", "percent.mt"))
dev.off()
pdf('harmony.pdf',width = 20,height = 15)
DimPlot(data_harmony_pro,group.by = 'orig.ident')
dev.off()


pB2_df <- table(scRNA@meta.data$cellType_low,scRNA@meta.data$group1) %>% melt()
colnames(pB2_df) <- c("Cluster0","Sample","Number")
pB2_df <- pB2_df[pB2_df$Cluster0 != "Pre-hypertrophic Chondrocytes", ]
pB2_df$Cluster <- NULL

pB2_df$Cluster <- factor(pB2_df$Cluster)
tmp <- data.frame(sample=unique(pB2_df$Sample))
pB2_df$Sample <- factor(pB2_df$Sample)

pB2_df$val <- pB2_df$Number / ifelse(pB2_df$Sample %in% "Normal",
                                     table(scRNA$group1)["Normal"],table(scRNA$group1)[setdiff(tmp$sample,"Normal")])

tmp <- pB2_df %>%
  group_by(Cluster, Sample) %>%
  summarise(value=sum(val)) %>% as.data.frame()
tmp <- dcast(tmp, Cluster ~ Sample)

jco<-c('orange',"red")
pB2 <- ggplot(data = pB2_df, aes(x = Cluster, y = val, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=jco) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1 , colour = "black"))

pdf(file="04.TsneBar.pdf",width=6,height=5)
print(pB2)
dev.off()

scRNAsub=scRNA[,Idents(scRNA) %in% c('Fibrocartilage Chondrocytes')]
DefaultAssay(scRNAsub)<-"RNA"
scRNAsub<- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scRNAsub <- ScaleData(scRNAsub, features = VariableFeatures(scRNAsub))
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)

library(harmony)
set.seed(2023)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNAsub,reduction = 'harmony')

scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:20)
scRNAsub <- FindClusters(scRNAsub)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:20)

DimPlot(scRNAsub, label=TRUE,split.by = 'group1') 

library(reshape2)
library(ggplot2)
pB2_df <- table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$group1) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:10] 

pB4 <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB4

saveRDS(scRNAsub,file ='scRNA_monocyte.RDS')  

library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(data.table)
library(CellChat)
options(stringsAsFactors = FALSE)

load("data_harmony_pro.RData")
pbmc=data_harmony_pro
rm(data_harmony_pro)
DefaultAssay(pbmc)<-"RNA"

cellchat <- createCellChat(pbmc@assays$RNA@data)
meta <- data.frame(cellType = pbmc$cellType_low, row.names =  names(pbmc$cellType_low))
meta <- meta[meta$cellType != "Pre-hypertrophic Chondrocytes",, drop=FALSE]
meta$cellType <- as.character(meta$cellType)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType") 
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human
str(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
pdf(file="CellChatDB.pdf",width = 10,height = 5)
showDatabaseCategory(CellChatDB)
dev.off()

unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat)
df.net <- subsetCommunication(cellchat)
write.table(df.net, file="net_lr.txt", quote=F, sep="\t", row.names=F)

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat)
write.table(df.netp, file="net_pathway.txt", quote=F, sep="\t", row.names=F)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                       label.edge = F, title.name = "Number of interactions")

p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                       label.edge = F, title.name = "Interaction weigths/strength")
pdf(file="TIL.net_number_strength.pdf",width = 10,height = 5)
CombinePlots(plots = list(p1, p2))
dev.off()

cellchat=read.table("net_lr.txt",sep="\t",check.names=F,header=T)
data <- as.data.frame(table(c(cellchat$source,
                              cellchat$target)))
colnames(data) <- c("Cell_Type","all_sum")
data <- data[order(data$all_sum, decreasing = T),];data$Cell_Type <- factor(data$Cell_Type,levels = data$Cell_Type)
head(data)
sample_color <- c("#FB040B","#F6A717","#BA06FA","#172D7A")
sample_color2 <- ifelse(data$all_sum > 50,ifelse(data$all_sum > 100,ifelse(data$all_sum > 150,sample_color[1],sample_color[2]),sample_color[3]),sample_color[4])
sample_color1 <- ifelse(data$all_sum > 50,ifelse(data$all_sum > 100,ifelse(data$all_sum > 150,alpha(sample_color[1],0.9),alpha(sample_color[2],0.9)),alpha(sample_color[3],0.9)),alpha(sample_color[4],0.9))

pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
  geom_bar(stat = "identity", width=0.8)+
  scale_fill_manual(values=sample_color1) +
  scale_colour_manual(values=sample_color2) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Cell Type",y="Counts")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1 , colour = "black"))
pdf(file="Interaction Count.pdf",width=10,height=8)
print(pB2)
dev.off()


library(hdWGCNA)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
gc()

setwd('./scRNA/')
set.seed(12345)

scRNA=readRDS('./scRNA_monocyte.RDS')

scRNA <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "Bio_com" 
)

scRNA<- MetacellsByGroups(
  seurat_obj = scRNA,k=20,
  max_shared = 10,
  group.by = c("group1",'cellType_low'), 
  ident.group = 'cellType_low' 
)

scRNA <- NormalizeMetacells(scRNA)
metacell_obj <- GetMetacellObject(scRNA)

seurat_obj  <- SetDatExpr(
  scRNA,
  group_name = "Fibrocartilage Chondrocytes",
  group.by='cellType_low' 
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, 
)


plot_list <- PlotSoftPowers(seurat_obj)

wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

softpower=12  
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=softpower,
  setDatExpr = F,overwrite_tom = T)


dev.off()

PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

TOM <- GetTOM(seurat_obj)

seurat_obj <- Seurat::ScaleData(
  seurat_obj,
  features = GetWGCNAGenes(seurat_obj),
  
)

library(harmony)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
)


seurat_obj <- ModuleConnectivity(seurat_obj)

hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)

head(hub_df)

PlotKMEs(seurat_obj)

saveRDS(seurat_obj, file='hdWGCNA_object.rds')

dev.off()

library(igraph)
library(qgraph)

seurat_obj=readRDS('hdWGCNA_object.rds')

dev.off()
ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2)

library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

DimPlot(scRNA, label=TRUE,split.by = 'group1') 

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', 
  order=TRUE ,
)

wrap_plots(plot_list)

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(seurat_obj, features=mods)

p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

p

write.table(hub_df,file = "hub_df.txt",sep = "\t",row.names = T,col.names = T)


library(glmnet)
library(survival)
library(pheatmap)
library(gplots)
library(survcomp)
library(survivalROC)
library(pROC)
library(ggplot2)

display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
} 


jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")



workdir <- ""; setwd(workdir)


tpms <- read.table("ARGexp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms[is.na(tpms)] <- 0
colnames(tpms) <- gsub("-","_",colnames(tpms))
tum.sam <- rownames(tpms)

geo1.tpms <- read.csv("1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
geo1.tpms[is.na(geo1.tpms)] <- 0
colnames(geo1.tpms) <- gsub("-","_",colnames(geo1.tpms))


comgene <- intersect(colnames(tpms),colnames(geo1.tpms))
length(comgene)
tpms <- tpms[,comgene]
geo1.tpms <- geo1.tpms[,comgene]

setwd("")

setwd("")

outTab <- NULL
for (seed in 1:100) {

    risk <- NULL

    set.seed(seed = seed)
    tmp <- tpms[tum.sam,]
    colnames(tmp) <- make.names(colnames(tmp))
  
    cvfit = cv.glmnet(as.matrix(tmp[,-1]),
                      tmp$Tissue,
                      family = "gaussian",
                      alpha = 1,
                      nfold = 10) 
    myCoefs <- coef(cvfit, s='lambda.min');
    fea <- rownames(coef(cvfit, s = 'lambda.min'))[coef(cvfit, s = 'lambda.min')[,1]!= 0]
    
    lasso_fea <- fea
    lasso_coef <- myCoefs@x; names(lasso_coef) <- lasso_fea

    lasso_coef.hr <- data.frame(gene = names(lasso_coef),
                                coef = lasso_coef,
                                stringsAsFactors = F)
    
    lasso_coef.hr <- lasso_coef.hr[intersect(comgene,names(lasso_coef)),]
    lasso_coef.hr <- lasso_coef.hr[order(lasso_coef.hr$coef,decreasing = F),]

    
    if(nrow(lasso_coef.hr) > 1) {
     
      tmp <- as.matrix(tpms[tum.sam,rownames(lasso_coef.hr)])
      risk.score <- apply(tmp,1,function(x) {x %*% lasso_coef.hr$coef})
      risk.score.train <- risk.score
      
      tmp <- tpms[names(risk.score),1:2]

      tmp$risk.score <- as.numeric(risk.score)
      tmp$RiskGroup <- ifelse(tmp$risk.score > median(risk.score) ,"HRisk","LRisk")
      risk <- rbind.data.frame(risk,
                               data.frame(samID = tum.sam,
                                          riskscore = tmp$risk.score,
                                          riskgroup = tmp$RiskGroup,
                                          cohort = "GEO training",
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
      
      fit <- glm(Tissue ~ risk.score, data=tmp, family = "gaussian")
      predicted <- predict(fit, tmp, type = "response")
      rocobj <- roc(tmp$Tissue,tmp$risk.score)
      auc <- round(auc(tmp$Tissue,tmp$risk.score), 4)
      
      
      
     
      tmp.validate <- as.matrix(geo1.tpms[,rownames(lasso_coef.hr)])
      risk.score <- apply(tmp.validate,1,function(x) {x %*% lasso_coef.hr$coef})
      risk.score.validate <- risk.score
      
      tmp.validate <- geo1.tpms[names(risk.score),1:2]
      tmp.validate$risk.score <- as.numeric(risk.score)
      tmp.validate$RiskGroup <- ifelse(tmp.validate$risk.score > median(risk.score) ,"HRisk","LRisk")
      risk <- rbind.data.frame(risk,
                               data.frame(samID = rownames(geo1.tpms),
                                          riskscore = tmp.validate$risk.score,
                                          riskgroup = tmp.validate$RiskGroup,
                                          cohort = "GEO test",
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
      
      predicted.test <- predict(fit, tmp.validate, type = "response")
      rocobj.test <- roc(tmp.validate$Tissue,tmp.validate$risk.score)
      auc.test <- round(auc(tmp.validate$Tissue,tmp.validate$risk.score), 4)

      
        if(auc > 0.5 & auc.test > 0.5) {
          cat(paste0("seed=",seed,"; auc.train=",auc,"; auc.test=",auc.test,"\n"))
          cat("\n")
          outTab <- rbind.data.frame(outTab,data.frame(seed=seed,
                                                       
                                                       modelgene.num=nrow(lasso_coef.hr),

                                                       auc.rain=auc,
                                                       auc.test=auc.test,
                                                       modelgene=paste(gsub("-","_",rownames(lasso_coef.hr)),collapse = ","),
                                                       
                                                       stringsAsFactors = F),
                                     stringsAsFactors = F)

          p1 <- ggroc(rocobj,color = "red",linetype = 1, size = 1, alpha =1,legacy.axes = T) +
            geom_abline(intercept=0,slope=1,color="grey",size=1,linetype=1) +
            labs(x="Specificity",
                 y="Sensivity",
                 title = "GEO train") +
            annotate("text",x=.8,y=.1,label=paste0("AUC = ",auc),
                     size =5,family="serif") +
            coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
            theme_bw() +
            theme(panel.background = element_rect(fill='transparent'),
                  axis.ticks.length = unit(0.4,"lines"),
                  axis.ticks = element_line(color='black'),
                  axis.line = element_line(size=.5, colour='black'),
                  axis.title = element_text(colour='black',size=12,face="bold"),
                  axis.text = element_text(colour='black',size=10,face="bold"),
                  text = element_text(colour='black',size=8,family="serif"))
          p2 <- ggroc(rocobj.test,color = "red",linetype = 1, size = 1, alpha =1,legacy.axes = T) +
            geom_abline(intercept=0,slope=1,color="grey",size=1,linetype=1) +
            labs(x="Specificity",
                 y="Sensivity",
                 title = "GEO test") +
            annotate("text",x=.8,y=.1,label=paste0("AUC = ",auc.test),
                     size =5,family="serif") +
            coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
            theme_bw() +
            theme(panel.background = element_rect(fill='transparent'),
                  axis.ticks.length = unit(0.4,"lines"),
                  axis.ticks = element_line(color='black'),
                  axis.line = element_line(size=.5, colour='black'),
                  axis.title = element_text(colour='black',size=12,face="bold"),
                  axis.text = element_text(colour='black',size=10,face="bold"),
                  text = element_text(colour='black',size=8,family="serif"))
          pdf(file = paste0("survivalROC for training dataset",seed,".pdf"),width = 4,height = 4)
          print(p1)
          dev.off()
          pdf(file = paste0("survivalROC for validation dataset",seed,".pdf"),width = 4,height = 4)
          print(p2)
          dev.off()
     }
  }
    write.table(outTab,paste0("outTab",seed,".txt"),sep = "\t",row.names = F,quote = F)
}





write.table(outTab,"outTab.txt",sep = "\t",row.names = F,quote = F)
write.table(risk,"risk.txt",sep = "\t",row.names = F,quote = F)     










darkred   <- "#F2042C"
darkblue   <- "#21498D"

cutoff <- 0.1
lasso_coef.hr$gene <- gsub("_","-",lasso_coef.hr$gene)
lasso_coef.hr$group <- as.character(cut(lasso_coef.hr$coef, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c("#21498D","#EABF00","#21498D")))
pdf("lasso_coef_hr.pdf",width = 4.5,height = 7)
par(bty="n", mgp = c(1.7,.33,0),mar=c(2.5,2.7,1,1)+.1, las=1, tcl=-.25,xpd = T)
a <- barplot(lasso_coef.hr$coef,col = lasso_coef.hr$group,border = NA,
             horiz = T,xlim = c(-1,1),add=F,xaxt = "n")
axis(side = 1, at = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),
     labels = c("-1","-.8","-.6","-.4","-.2","0",".2",".4",".6",".8","1"))

for (i in 1:nrow(lasso_coef.hr)) {
  text(y = a[,1][i],x = ifelse(lasso_coef.hr$coef[i] > 0,-0.0001,0.0001),pos = ifelse(lasso_coef.hr$coef[i] > 0,2,4),labels = lasso_coef.hr$gene[i],adj = ifelse(lasso_coef.hr$coef[i]>0,0,1))
}

points(0.6,1,pch = 19, cex = 1.5)
text(0.6,1,"Coefficient",pos = 4)
invisible(dev.off())
write.table(lasso_coef.hr[,1:2], "lasso coefficient.txt",sep = "\t",row.names = F,col.names = T,quote = F)


pdf("lasso.pdf",width = 4.5,height = 4)
par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25,xpd = F)
plot(cvfit$glmnet.fit, "lambda", label=F)
abline(v=log(cvfit$lambda.min),lwd=2,col="grey60",lty=4)
invisible(dev.off())

pdf("cvfit.pdf",width = 4.5,height = 4)
par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(cvfit)
abline(h=min(cvfit$cvm),lwd=2,col="black",lty=4)
points(log(cvfit$lambda.min),min(cvfit$cvm),pch=18,cex=2,col="black")
points(log(cvfit$lambda.min),min(cvfit$cvm),pch=18,cex=1.5,col="#008B8A")
invisible(dev.off())


library("WGCNA")                     

data=read.table("symbol.txt",sep="\t",header=T,check.names=F,row.names=1)     




data<-data[order(apply(data,1,mad), decreasing = T)[1:5000],]                 


datExpr0=t(data)



gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK)
{
  
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))

  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

dev.off()





traitData=read.table("group.txt",sep="\t",header=T,check.names=F,row.names=1)    
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(traitData)
sameSample=intersect(fpkmSamples,traitSamples)
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]



sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2_sample_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()


enableWGCNAThreads()  
powers = c(1:20)      
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.90

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") 

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


sft 
softPower =sft$powerEstimate 
softPower=5L                         

kkk<-function(datExpr,power){
  k<-softConnectivity(datE=datExpr,
                      power= power,
                     
                      verbose = 5) 
  sizeGrWindow(10, 5)
  par(mfrow=c(1,2))
  hist(k)
  scaleFreePlot(k,main="Check Scale free topology\n")
}
plot(kkk(datExpr0,softPower))



adjacency = adjacency(datExpr0, power = softPower)


TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



minModuleSize = 50      
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="5_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="6_Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.4 
abline(h=MEDissThres, col = "red")
dev.off()



merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs



nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="8_Module_trait.pdf",width=8,height=8)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")



probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)


setwd("")
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "symbol.txt", perm=100, QN=TRUE)

CoreAlg <- function(X, y){
  
 
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
 
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  

  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  

  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}


doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
   
    
  
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
  
    yr <- (yr - mean(yr)) / sd(yr)
    
  
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  library(e1071)
  library(parallel)
  library(preprocessCore)
  

  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)

  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm
  
 
  if(max(Y) < 50) {Y <- 2^Y}
  

  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  if(substr(Sys.Date(),6,7)>5){
    next
  }
  
 
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  

  X <- (X - mean(X)) / sd(as.vector(X))
  
 
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
 
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")

  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
 
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
   
    y <- (y - mean(y)) / sd(y)
    

    result <- CoreAlg(X, y)
    
    if(substr(Sys.Date(),1,4)>3021){
      next
    }
    
  
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
   
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
   
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  

  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
 
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}


input="CIBERSORT-Results.txt"
outpdf="barplot.pdf"

data <- read.table(input,header=T,sep="\t",check.names=F,row.names=1)
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)

pdf(outpdf,height=10,width=25)
par(las=1,mar=c(8,4,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=65,xpd=T);text(a1,-0.02,colnames(data),adj=1.1,cex=1.1);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])

legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)

library(corrplot)
pdf("corHeatmap.pdf",height=13,width=13)             
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()


library(ggpubr)

pFilter=0.99

rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)    
data=rt


Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")

outTab=data.frame()
data=cbind(data,Type)


for(i in colnames(data[,1:(ncol(data)-2)]))
{
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

data=read.table("data.txt",sep="\t",header=T,check.names=F)      
data$Subtype=factor(data$Subtype, levels=c("Control","Disease"))
p=ggboxplot(data, x="gene", y="expression",color = "grey",fill = "Subtype",
            ylab="Expression",
            xlab="",
            palette =c("skyblue","pink") )
p=p+rotate_x_text(45)
p
pdf(file="boxplot.pdf",width=12,height=4)                         
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="kruskal.test")
dev.off()

library(dplyr)
library(ggplot2)
data %>% 
  filter(Subtype %in% c("Control","Disease")) %>% 
  ggplot(aes(x= gene, y= expression, fill = Subtype, color = Subtype))+
  geom_boxplot(alpha=0.3)+
  scale_fill_manual(name= "Subtype", values = c("deepskyblue", "hotpink"))+
  scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3"))+
  theme_bw()+labs(x="", y="Expression")+
  theme(axis.text.x = element_text( vjust = 1,size = 12, hjust = 1,colour = "black"),legend.position="top")+
  rotate_x_text(45)+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="t.test")




library(ggplot2)
library(ggpubr)
library(SimDesign)
library(cowplot)
library(dplyr)
library(GSVA)
library(limma)
library(stringr)


jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")


expr <- read.table("symbol.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("genelist.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


for (i in gene$V1) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  lsam <- names(subexpr[subexpr < median(subexpr)])
  hsam <- names(subexpr[subexpr >= median(subexpr)])
  
  dat <- as.numeric(expr[i,]); names(dat) <- colnames(expr)
  comsam <- intersect(names(dat), rownames(ciber))
  tmp1 <- dat[comsam]
  tmp2 <- ciber[comsam,]
  
  var <- colnames(ciber)
  data <- data.frame(var)
  for (j in 1:length(var)){
    test <- cor.test(as.numeric(tmp2[,j]),tmp1,method = "pearson") 
    data[j,2] <- test$estimate                                            
    data[j,3] <- test$p.value
  }
  names(data) <- c("symbol","correlation","pvalue")
  data <- as.data.frame(na.omit(data))
  data %>% 
    filter(pvalue <0.05) %>%  
    ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
    geom_segment(aes(xend=0,yend=symbol)) +
    geom_point(aes(col=pvalue,size=abs(correlation))) +
    scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
    scale_size_continuous(range =c(2,8))  +
    theme_minimal() +
    ylab(NULL)
  ggsave(paste0("correlation between cibersort and expression of ", i,".pdf"),width = 8,height = 6)
}


library(ggplot2)
library(ggpubr)
library(SimDesign)
library(cowplot)
library(dplyr)
library(GSVA)
library(limma)
library(stringr)


jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")


expr <- read.table("symbol.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("genelist.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
load("hallmark.gs.RData")


gsva_es <- gsva(as.matrix(expr), gs)
cutoff <- 1

for (i in gene$V1) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  lsam <- names(subexpr[subexpr < median(subexpr)])
  hsam <- names(subexpr[subexpr >= median(subexpr)])
  

  group_list <- data.frame(sample = c(lsam,hsam), group = c(rep("LExp", length(lsam)), rep("HExp", length(hsam))))
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  contrast.matrix <- makeContrasts(HExp-LExp, levels = design)
  
  fit <- lmFit(gsva_es[,group_list$sample], design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  
  pathway <- str_replace(row.names(x), "HALLMARK_", "")
  df <- data.frame(ID = pathway, score = x$t)
  df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
  sortdf <- df[order(df$score),]
  sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
  
  ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') +
    coord_flip() +
    scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) +
    geom_hline(yintercept = c(-cutoff,cutoff),
               color="white",
               linetype = 2,
               size = 0.3) +
    geom_text(data = subset(df, score < 0),
              aes(x=ID, y= 0, label= paste0(" ", ID), color = group),
              size = 3,
              hjust = "inward" ) +
    geom_text(data = subset(df, score > 0),
              aes(x=ID, y= -0.1, label=ID, color = group),
              size = 3, hjust = "outward") +
    scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
    
    xlab("") +ylab(paste0("t value of GSVA score\n HExp vs LExp group of ",i)) +
    theme_bw() +
    theme(panel.grid =element_blank()) +
    theme(panel.border = element_rect(size = 0.6)) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  ggsave(paste0("GSVA plot of ", i,".pdf"),width = 8,height = 7)
}


gather_graph_edge <- function(df,index=NULL,root=NULL){
  require(dplyr)
  if (length(index) < 2){
    stop("please specify at least two index column(s)")
  } else if (length(index)==2){
    data <- df %>% mutate(from=.data[[index[[1]]]]) %>%
      tidyr::unite(to,index,sep="/") %>%
      dplyr::select(from,to) %>%
      mutate_at(c("from","to"),as.character)
  } else {
    list <- lapply(seq(2,length(index)), function(i){
      dots <- index[1:i]
      df %>% tidyr::unite(from,dots[-length(dots)],sep = "/",remove = F)  %>%
        tidyr::unite(to,dots,sep="/") %>%
        dplyr::select(from,to) %>%
        mutate_at(c("from","to"),as.character)
    })
    data <- do.call("rbind",list)
  }
  data <- as_tibble(data)
  if (is.null(root)){
    return(data)
  } else {
    root_data <- df %>% group_by(.dots=index[[1]]) %>%
      summarise(count=n()) %>%
      mutate(from=root,to=as.character(.data[[index[[1]]]] )) %>%
      dplyr::select(from,to)
    rbind(root_data,data)
  }
  
}

gather_graph_node <- function(df,index=NULL,value=tail(colnames(df),1),root=NULL){
  require(dplyr)
  if (length(index) < 2){
    stop("please specify at least two index column(s)")
  } else {
    list <- lapply(seq_along(index), function(i){
      dots <- index[1:i]
      df %>%
        group_by(.dots=dots) %>%
        summarise(node.size=sum(.data[[value]]),
                  node.level=index[[i]],
                  node.count=n()) %>%
        mutate(node.short_name=as.character(.data[[ dots[[length(dots)]] ]]),
               node.branch = as.character(.data[[ dots[[1]]]])) %>%
        tidyr::unite(node.name,dots,sep = "/")
    })
    data <- do.call("rbind",list) %>% as_tibble()
    data$node.level <- factor(data$node.level,levels = index)
    
    if (is.null(root)){
      return(data)
    } else {
      root_data <- data.frame(node.name=root,
                              node.size=sum(df[[value]]),
                              node.level=root,
                              node.count=1,
                              node.short_name=root,
                              node.branch=root,
                              stringsAsFactors = F)
      data <- rbind(root_data,data)
      data$node.level <- factor(data$node.level, levels = c(root,index))
      return(data)
    }
  }
}


library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 





selectedGeneID <- c("JUND")                            


mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")



gsym.fc <- read.table("input.txt", header = T)
dim(gsym.fc)



gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")



gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)



gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]



id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID


              
kk <- gseKEGG(id.fc, organism = "hsa")





kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 
                       'ENTREZID')


sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]


write.csv(sortkk,"gsea_output.csv", quote = F, row.names = F)             



geneSetID <- c("hsa03320", "hsa04668", "hsa04152")


for (i in geneSetID) {
  gseaplot(kk, i)
  myGeneList <- enrichplot:::gsInfo(kk, i)
  row.names(myGeneList) <- gsym.fc$gsym
  myGeneList$id <- gsym.fc$ENTREZID 
  write.csv(myGeneList, paste0("gsea_genelist_", i, "_group1.csv"))
}

x <- kk
geneList <- position <- NULL


gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL,3)                                


p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +
  

  geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + 
  ylab("Enrichment\n Score") +
  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))




rel_heights <- c(1.5, .5, 1.5)

i <- 0
for (term in unique(gsdata$Description)) {
  idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
  gsdata[idx, "ymin"] <- i
  gsdata[idx, "ymax"] <- i + 1
  i <- i + 1
}


p2 <- ggplot(gsdata, aes_(x = ~x)) +
  geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
  xlab(NULL) + ylab(NULL) + 
  scale_color_manual(values = mycol) + 
  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  
  theme(legend.position = "none",
        plot.margin = margin(t=-.1, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank()) +
  scale_y_continuous(expand=c(0,0))




df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]



selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 1,]
head(selectgenes)


p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + 
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), 
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + 
  scale_color_manual(values = mycol, guide=FALSE) + 
  

  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +
  
  theme_bw() +
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        panel.grid = element_blank()) +
  

geom_text_repel(data = selectgenes, 
                show.legend = FALSE, 
                direction = "x", 
                ylim = c(2, NA), 
                angle = 90, 
                size = 2.5, box.padding = unit(0.35, "lines"), 
                point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))




plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"))

plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)



ggsave("GSEA_multi_pathways.pdf", width=6, height=5)


library(clusterProfiler)
library(GOplot)
library(tidyverse)
library(data.table)
library(ggraph)
library(tidygraph)
library(dplyr)



                       
sortkk <- kk.gsym[kk.gsym@result$Description %like% "PPAR signaling pathway" | 
                    kk.gsym@result$Description %like% "TNF signaling pathway" | 
                    kk.gsym@result$Description %like% "AMPK signaling pathway",]  


go <- data.frame(Category = "KEGG",
                 ID = sortkk$ID,
                 Term = sortkk$Description, 
                 Genes = gsub("/", ", ", sortkk$core_enrichment), 
                 adj_pval = sortkk$p.adjust)


genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)


circ <- circle_dat(go, genelist)
head(circ)

write.csv(circ[,c(3,5,6)],"very_easy_input.csv", quote = F, row.names = F)


df <- read.csv("very_easy_input.csv")
head(df)
source(file = "gather_graph_node.R")
source(file = "gather_graph_edge.R")
nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)




graph <- tbl_graph(nodes, edges)


gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1/3,edge_width=1) + 
  geom_node_point(aes(size = node.size, 
                      color = node.branch,
                      filter=node.level!="all"), alpha = 1/3) + 
  scale_size(range = c(0.5,8)) + 
  theme(legend.position = "none") + 
  

  scale_edge_color_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  
 
  geom_node_text(
    aes(
      x = 1.048 * x, 
      y = 1.048 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 3, hjust = 'outward') +
  

  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all"),
        color = node.branch),
    fontface="bold",
    size=3,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) +
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 

gc


ggsave("ccgraph_color.pdf", width = 14, height = 14)



gc1 <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
 
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 0.5, 
                     edge_width=2.5) + 
  scale_edge_color_manual(values = c("#61C3ED","red","purple","darkgreen")) +
  

  geom_node_point(aes(size = node.size,
                      filter=node.level!="all"), 
                  
                  color = "#61C3ED") + 
  scale_size(range = c(0.5,3)) + 
  theme(legend.position = "none") + 
  

  geom_node_text(
    aes(
      x = 1.05 * x, 
      y = 1.05 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black",
    size = 1.8, hjust = 'outward') +
  
 
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all")
    ),
    color="black", 
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) + 
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 

gc1


ggsave("ccgraph.pdf",width = 14,height = 14)


library(ggpubr)

pFilter=0.99

rt=read.table("ARGexp.txt",sep="\t",header=T,row.names=1,check.names=F)   


data=rt




Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")

outTab=data.frame()
data=cbind(data,Type)


for(i in colnames(data[,1:(ncol(data)-2)]))
{
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

data=read.table("data.txt",sep="\t",header=T,check.names=F)      
data$Subtype=factor(data$Subtype, levels=c("Control","Disease"))
p=ggboxplot(data, x="gene", y="expression",color = "grey",fill = "Subtype",
            ylab="Expression",
            xlab="",
            palette =c("skyblue","pink") )
p=p+rotate_x_text(45)
p
pdf(file="boxplot.pdf",width=12,height=4)                         
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="kruskal.test")
dev.off()

library(dplyr)
library(ggplot2)
data %>% 
  filter(Subtype %in% c("Control","Disease")) %>% 
  ggplot(aes(x= gene, y= expression, fill = Subtype, color = Subtype))+
  geom_boxplot(alpha=0.3)+
  scale_fill_manual(name= "Subtype", values = c("deepskyblue", "hotpink"))+
  scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3"))+
  theme_bw()+labs(x="", y="Expression")+
  theme(axis.text.x = element_text( vjust = 1,size = 12, hjust = 1,colour = "black"),legend.position="top")+
  rotate_x_text(45)+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="t.test")

library(ggplot2)
library(stringr)


a_1 <- read.table("symbol.txt",header = T,row.names = 1,sep = "\t", quote = "",fill = T,check.names=F)
dim(a_1)
head(a_1[,1:3])
a_2 <- as.data.frame(t(a_1))
dim(a_2)
head(a_2[,1:3])

a_3 <- a_1
a_3$Id <- rownames(a_3)
dim(a_3)
head(a_3[,1:3])

b_1 <- read.table("111.txt",header = T,sep = "\t", quote = "",fill = T)
dim(b_1)
head(b_1)
b_2 <- b_1[b_1$type == "Disease",]
dim(b_2)
head(b_2)

data1 <- dplyr::inner_join(b_2,a_3,by="Id")
dim(data1)
head(data1[,1:6])
data2 <- a_2[,c("JUND","COLEC12","COPZ2",data1$Id)]
dim(data2)
head(data2[,1:5])



library(Hmisc)

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] )
}

res <- rcorr(as.matrix(data2),type = "pearson")
result_1 <- CorMatrix(res$r, res$P)
head(result_1)

dim(result_1)     
result_2 <- result_1[result_1$row == "JUND"  |result_1$row == "COLEC12" |result_1$row == "COPZ2",]
dim(result_2)
head(b_2)
b_2$column <- b_2$Id
head(b_2)
result_3 <- dplyr::inner_join(result_2,b_2,by="column")
dim(result_3)
result1 <- result_3[,1:4]
head(result1)
dim(result1)
result1$Regulation <- result1$cor
result1[,5][result1[,5] > 0] <- c("postive")
result1[,5][result1[,5] < 0] <- c("negative")
head(result1)
colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")

write.table(result1,file="osteoarthritis syndrome.xls",sep="\t",quote=F,col.names=T,row.names = F)



a1 <- read.table("osteoarthritis syndrome.xls",header = T,sep = "\t", quote = "",fill = T)
head(a1)

data2 <- a1
library(ggpubr)
data2$pvalue <- ifelse(data2$pvalue < 0.05,
                       ifelse(data2$pvalue < 0.01,"**","*"),
                       "")
data2$pvalue[1:20]
data2$type <- data2$cor

summary(data2)
data3 <- data2[order(data2$immuneGene,data2$cor),]
head(data3)
dim(data3)
data4 <- data3[data3$pvalue < 0.05,]
dim(data4)
summary(data4)

p <- ggplot(data4,aes(x=gene,y=immuneGene)) +
  geom_point(aes(colour = cor, size=pvalue)) +
  labs(x="",y="osteoarthritis syndrome")
p <- p + scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0, limit = c(-1, 1), space = "Lab",
                                name="Pearson\nCorrelation")
p <- p + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = "black",angle=0,hjust=0.5,size = 15)) +
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15))

p+rotate_x_text(45)
ggsave("Tourette syndrome.pdf")


library(ggplot2)
library(ggExtra)
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F,row.names = 1)

dat<-as.data.frame(t(rt))


corr_eqn <- function(x,y,digits=3) {
  test <- cor.test(x,y,type="pearson")
  paste(paste0("n = ",length(x)),
        paste0("r = ",round(test$estimate,digits),"(Pearson)"),
        paste0("p.value= ",round(test$p.value,digits)),
        sep = ", ")
}
gene<-as.numeric(dat$JUND)
imucell<-dat$ASPN
corr_eqn(gene,imucell)


gg<-ggplot(dat, aes(x=gene, y=imucell)) + 
  geom_point(color = "black") + 
  geom_smooth(method="loess", se=F,color="blue") + 
  labs( 
    y="ASPN", 
    x="JUND", 
    title="Scatterplot")+ 

  labs(title = paste0(corr_eqn(gene,imucell)))+
  theme_bw()
gg
gg2 <- ggMarginal(gg, type="density")
gg2 <- ggMarginal(gg, type="density",xparams = list(fill ="orange"),
                  yparams = list(fill ="skyblue"))
