library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)

NT.data <- Read10X(data.dir ="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/1_NT_filtered_feature_bc_matrix")
NT <- CreateSeuratObject(counts =  NT.data$`Gene Expression`, project = "NT")
NT [["percent.mt"]] <- PercentageFeatureSet(object =  NT, pattern = "^mt-")
VlnPlot(object =  NT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
plot1 <- FeatureScatter(object =  NT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object =  NT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#subset data, for now only exclude cells with %mito above 10%
#NT <- subset(x =  NT, subset = percent.mt < 10)
# normalize the data
NT <- NormalizeData(object =  NT, normalization.method = "LogNormalize", scale.factor = 10000)
NT <- FindVariableFeatures(object =  NT, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  NT)
NT <- ScaleData(object = NT, features = all.genes)
#NT <- SCTransform(NT, vars.to.regress = "percent.mt", verbose = FALSE)
NT <- RunPCA(object =  NT, features = VariableFeatures(object =  NT))
ElbowPlot(object =  NT, ndims = 50)
NT <- FindNeighbors(object =  NT, dims = 1:20)
NT <- FindClusters(object =  NT, resolution = 1.0)
#proceed to UMAP:
NT <- RunUMAP(object =  NT, dims = 1:20)
DimPlot(object =  NT, reduction = "umap", label = TRUE)
FeaturePlot(object= NT, features = c("Slc17a6"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')

#subset just neurons
NT_NeuronsOnly <- subset(x = NT, idents = c(4, 11, 0, 8), invert = TRUE)
#NT_NeuronsOnly <- SCTransform(NT_NeuronsOnly, vst.flavor = "v2", verbose = FALSE)
NT_NeuronsOnly <- FindVariableFeatures(object =  NT_NeuronsOnly, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  NT_NeuronsOnly)
NT_NeuronsOnly <- ScaleData(object =  NT_NeuronsOnly, features = all.genes)
NT_NeuronsOnly <- RunPCA(object =  NT_NeuronsOnly, features = VariableFeatures(object =  NT_NeuronsOnly))
ElbowPlot(object =  NT_NeuronsOnly, ndims = 50)
NT_NeuronsOnly <- FindNeighbors(object =  NT_NeuronsOnly, dims = 1:20)
NT_NeuronsOnly <- FindClusters(object =  NT_NeuronsOnly, resolution = 1.0) #2.0 before
NT_NeuronsOnly <- RunUMAP(object =  NT_NeuronsOnly, dims = 1:20, return.model=TRUE)
DimPlot(object =  NT_NeuronsOnly, reduction = "umap", label=TRUE, pt.size = 0.25)
#FeaturePlot(object= NT_NeuronsOnly, features = c("Npy1r"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')

save(NT_NeuronsOnly, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/1_NT/NT.RData")
rm(NT.data, NT)

## Repeat this for 2_PBS conditions:

PBS.data <- Read10X(data.dir ="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/2_PBS_filtered_feature_bc_matrix")
PBS <- CreateSeuratObject(counts =  PBS.data$`Gene Expression`, project = "PBS")
PBS [["percent.mt"]] <- PercentageFeatureSet(object =  PBS, pattern = "^mt-")
VlnPlot(object =  PBS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
plot1 <- FeatureScatter(object =  PBS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object =  PBS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# normalize the data
PBS <- NormalizeData(object =  PBS, normalization.method = "LogNormalize", scale.factor = 10000)
PBS <- FindVariableFeatures(object =  PBS, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  PBS)
PBS <- ScaleData(object =  PBS, features = all.genes)
PBS <- RunPCA(object =  PBS, features = VariableFeatures(object =  PBS))
ElbowPlot(object =  PBS, ndims = 50)
PBS <- FindNeighbors(object =  PBS, dims = 1:20)
PBS <- FindClusters(object =  PBS, resolution = 1.0)
#proceed to UMAP:
PBS <- RunUMAP(object =  PBS, dims = 1:20)
DimPlot(object =  PBS, reduction = "umap", label = TRUE)
FeaturePlot(object= PBS, features = c("Slc17a6"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')##

PBS_NeuronsOnly <- subset(x = PBS, idents = c(4,5,0,18), invert = TRUE)
#PBS_NeuronsOnly <- SCTransform(PBS_NeuronsOnly, vst.flavor = "v2", verbose = FALSE)
PBS_NeuronsOnly <- FindVariableFeatures(object =  PBS_NeuronsOnly, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  PBS_NeuronsOnly)
PBS_NeuronsOnly <- ScaleData(object =  PBS_NeuronsOnly, features = all.genes)
PBS_NeuronsOnly <- RunPCA(object =  PBS_NeuronsOnly, features = VariableFeatures(object =  PBS_NeuronsOnly))
ElbowPlot(object =  PBS_NeuronsOnly, ndims = 50)
PBS_NeuronsOnly <- FindNeighbors(object =  PBS_NeuronsOnly, dims = 1:20)
PBS_NeuronsOnly <- FindClusters(object =  PBS_NeuronsOnly, resolution = 1.0)
PBS_NeuronsOnly <- RunUMAP(object =  PBS_NeuronsOnly, dims = 1:20)
DimPlot(object =  PBS_NeuronsOnly, reduction = "umap", label=TRUE, pt.size = 0.25)
#FeaturePlot(object= PBS, features = c("P2ry1"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')##

save(PBS_NeuronsOnly, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/2_PBS/PBS.RData")
rm(PBS.data, PBS)

## Repeat this for 3_OVA conditions:

OVA.data <- Read10X(data.dir ="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/3_OVA_filtered_feature_bc_matrix")
OVA <- CreateSeuratObject(counts =  OVA.data$`Gene Expression`, project = "OVA")
OVA[["percent.mt"]] <- PercentageFeatureSet(object =  OVA, pattern = "^mt-")
VlnPlot(object =  OVA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
plot1 <- FeatureScatter(object =  OVA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object =  OVA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# normalize the data
OVA <- NormalizeData(object =  OVA, normalization.method = "LogNormalize", scale.factor = 10000)
OVA <- FindVariableFeatures(object =  OVA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  OVA)
OVA <- ScaleData(object =  OVA, features = all.genes)
OVA <- RunPCA(object =  OVA, features = VariableFeatures(object =  OVA))
ElbowPlot(object =  OVA, ndims = 50)
OVA <- FindNeighbors(object =  OVA, dims = 1:20)
OVA<- FindClusters(object =  OVA, resolution = 1.0)
#proceed to UMAP:
OVA <- RunUMAP(object =  OVA, dims = 1:20)
DimPlot(object =  OVA, reduction = "umap", label = TRUE)
FeaturePlot(object= OVA, features = c("Slc17a6"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')##

#Subset just neurons:
OVA_NeuronsOnly <- subset(x = OVA, idents = c(3,0,7), invert = TRUE)
#OVA_NeuronsOnly <- SCTransform(OVA_NeuronsOnly, vst.flavor = "v2", verbose = FALSE)
OVA_NeuronsOnly <- FindVariableFeatures(object =  OVA_NeuronsOnly, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  OVA_NeuronsOnly)
OVA_NeuronsOnly <- ScaleData(object =  OVA_NeuronsOnly, features = all.genes)
OVA_NeuronsOnly <- RunPCA(object =  OVA_NeuronsOnly, features = VariableFeatures(object =  OVA_NeuronsOnly))
ElbowPlot(object =  OVA_NeuronsOnly, ndims = 50)
OVA_NeuronsOnly <- FindNeighbors(object =  OVA_NeuronsOnly, dims = 1:20)
OVA_NeuronsOnly <- FindClusters(object =  OVA_NeuronsOnly, resolution = 1.0)
OVA_NeuronsOnly <- RunUMAP(object =  OVA_NeuronsOnly, dims = 1:20)
DimPlot(object =  OVA_NeuronsOnly, reduction = "umap", label=TRUE, pt.size = 0.25)

save(OVA_NeuronsOnly, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/3_OVA/OVA.RData")
rm(OVA.data, OVA)

## Repeat this for 4_PBS_Chase conditions:

PBSChase.data <- Read10X(data.dir ="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/4_PBS_C_filtered_feature_bc_matrix")
PBSChase <- CreateSeuratObject(counts =  PBSChase.data$`Gene Expression`, project = "PBSChase")
PBSChase[["percent.mt"]] <- PercentageFeatureSet(object =  PBSChase, pattern = "^mt-")
VlnPlot(object =  PBSChase, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
plot1 <- FeatureScatter(object =  PBSChase, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object =  PBSChase, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# normalize the data
PBSChase <- NormalizeData(object =  PBSChase, normalization.method = "LogNormalize", scale.factor = 10000)
PBSChase <- FindVariableFeatures(object =  PBSChase, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  PBSChase)
PBSChase <- ScaleData(object =  PBSChase, features = all.genes)
PBSChase <- RunPCA(object =  PBSChase, features = VariableFeatures(object =  PBSChase))
ElbowPlot(object =  PBSChase, ndims = 50)
PBSChase <- FindNeighbors(object =  PBSChase, dims = 1:20)
PBSChase<- FindClusters(object =  PBSChase, resolution = 1.0)
#proceed to UMAP:
PBSChase <- RunUMAP(object =  PBSChase, dims = 1:20)
DimPlot(object =  PBSChase, reduction = "umap", label = TRUE)
FeaturePlot(object= PBSChase, features = c("Slc17a6"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')##

#Subset just neurons:
PBSChase_NeuronsOnly <- subset(x = PBSChase, idents = c(1,10,5), invert = TRUE)
#PBSChase_NeuronsOnly <- SCTransform(PBSChase_NeuronsOnly, vst.flavor = "v2", verbose = FALSE)
PBSChase_NeuronsOnly <- FindVariableFeatures(object =  PBSChase_NeuronsOnly, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  PBSChase_NeuronsOnly)
PBSChase_NeuronsOnly <- ScaleData(object =  PBSChase_NeuronsOnly, features = all.genes)
PBSChase_NeuronsOnly <- RunPCA(object =  PBSChase_NeuronsOnly, features = VariableFeatures(object =  PBSChase_NeuronsOnly))
ElbowPlot(object =  PBSChase_NeuronsOnly, ndims = 50)
PBSChase_NeuronsOnly <- FindNeighbors(object =  PBSChase_NeuronsOnly, dims = 1:20)
PBSChase_NeuronsOnly <- FindClusters(object =  PBSChase_NeuronsOnly, resolution = 1.0)
PBSChase_NeuronsOnly <- RunUMAP(object =  PBSChase_NeuronsOnly, dims = 1:20)
DimPlot(object =  PBSChase_NeuronsOnly, reduction = "umap", label=TRUE, pt.size = 0.25)

save(PBSChase_NeuronsOnly, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/4_PBS_C/PBSChase.RData")
rm(PBSChase.data, PBSChase)

## Repeat this for 5_OVA_Chase conditions:

OVAChase.data <- Read10X(data.dir ="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/5_OVA_C_filtered_feature_bc_matrix")
OVAChase <- CreateSeuratObject(counts =  OVAChase.data$`Gene Expression`, project = "OVAChase")
OVAChase[["percent.mt"]] <- PercentageFeatureSet(object =  OVAChase, pattern = "^mt-")
VlnPlot(object =  OVAChase, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
plot1 <- FeatureScatter(object =  OVAChase, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object =  OVAChase, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# normalize the data
OVAChase <- NormalizeData(object =  OVAChase, normalization.method = "LogNormalize", scale.factor = 10000)
OVAChase <- FindVariableFeatures(object =  OVAChase, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  OVAChase)
OVAChase <- ScaleData(object =  OVAChase, features = all.genes)
OVAChase <- RunPCA(object =  OVAChase, features = VariableFeatures(object =  OVAChase))
ElbowPlot(object =  OVAChase, ndims = 50)
OVAChase <- FindNeighbors(object =  OVAChase, dims = 1:20)
OVAChase<- FindClusters(object =  OVAChase, resolution = 1.0)
#proceed to UMAP:
OVAChase <- RunUMAP(object =  OVAChase, dims = 1:20)
DimPlot(object =  OVAChase, reduction = "umap", label = TRUE)
FeaturePlot(object= OVAChase, features = c("Slc17a6"), cols = c("grey","blue"), order = T, min.cutoff = 0, max.cutoff = 'q99')##

#Subset just neurons:
OVAChase_NeuronsOnly <- subset(x = OVAChase, idents = c(1,3,5), invert = TRUE)
#OVAChase_NeuronsOnly <- SCTransform(OVAChase_NeuronsOnly, vst.flavor = "v2", verbose = FALSE)
OVAChase_NeuronsOnly <- FindVariableFeatures(object =  OVAChase_NeuronsOnly, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  OVAChase_NeuronsOnly)
OVAChase_NeuronsOnly <- ScaleData(object =  OVAChase_NeuronsOnly, features = all.genes)
OVAChase_NeuronsOnly <- RunPCA(object =  OVAChase_NeuronsOnly, features = VariableFeatures(object =  OVAChase_NeuronsOnly))
ElbowPlot(object =  OVAChase_NeuronsOnly, ndims = 50)
OVAChase_NeuronsOnly <- FindNeighbors(object =  OVAChase_NeuronsOnly, dims = 1:20)
OVAChase_NeuronsOnly <- FindClusters(object =  OVAChase_NeuronsOnly, resolution = 1.0)
OVAChase_NeuronsOnly <- RunUMAP(object =  OVAChase_NeuronsOnly, dims = 1:20, return.model = T)
DimPlot(object =  OVAChase_NeuronsOnly, reduction = "umap", label=TRUE, pt.size = 0.25)

save(OVAChase_NeuronsOnly, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/5_OVA_C/OVAChase.RData")
rm(OVAChase.data, OVAChase)


#Import Data

load("C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/1_NT/NT.RData")
load("C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/2_PBS/PBS.RData")
load("C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/3_OVA/OVA.RData")
load("C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/4_PBS_C/PBSChase.RData")
load("C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/5_OVA_C/OVAChase.RData")



##Next, integrate all the datasets
NT_NeuronsOnly$condition = "NT"
NT_NeuronsOnly$condition2 ="control"
PBS_NeuronsOnly$condition = "PBS"
PBS_NeuronsOnly$condition2 = "control"
OVA_NeuronsOnly$condition = "OVA"
OVA_NeuronsOnly$condition2 = "OVA-treated"
PBSChase_NeuronsOnly$condition = "PBS_Chase"
PBSChase_NeuronsOnly$condition2 = "control"
OVAChase_NeuronsOnly$condition = "OVA_Chase"
OVAChase_NeuronsOnly$condition2 = "OVA-treated"

integrate.list <- list(NT_NeuronsOnly, PBS_NeuronsOnly, OVA_NeuronsOnly, PBSChase_NeuronsOnly, OVAChase_NeuronsOnly)
new.anchors <- FindIntegrationAnchors(object.list = integrate.list, dims = 1:30, scale = T, normalization.method = "LogNormalize")
multiome.combined <- IntegrateData(anchorset = new.anchors, dims = 1:30,  normalization.method = "LogNormalize" )
#alternatively, use SCT -- but overregresses
#features <- SelectIntegrationFeatures(object.list = integrate.list, nfeatures = 3000)
#integrate.list <- PrepSCTIntegration(object.list = integrate.list, anchor.features = features)
#multiome.combined <- SCTransform(multiome.combined, vst.flavor = "v2", verbose = FALSE)
#multiome.combined <- FindVariableFeatures(object =  multiome.combined, selection.method = "vst", nfeatures = 2000)

multiome.combined[["RNA"]] <- JoinLayers(multiome.combined[["RNA"]])
DefaultAssay(object = multiome.combined) <- "RNA"
multiome.combined <- NormalizeData(object =  multiome.combined, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x =  multiome.combined)
multiome.combined <- ScaleData(object = multiome.combined, features = all.genes, model.use = "LogNormalize", verbose = FALSE, var.to.regress = "percent.mt")
multiome.combined <- RunPCA(object = multiome.combined, npcs = 50, verbose = FALSE)
ElbowPlot(object = multiome.combined, ndims = 50)
multiome.combined <- FindNeighbors(multiome.combined, dims = 1:30)
multiome.combined <- FindClusters(multiome.combined, resolution = 1.0)
multiome.combined <- RunUMAP(object = multiome.combined, reduction = "pca", dims = 1:30)
DimPlot(object =  multiome.combined, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.25)

#FeaturePlot(object= multiome.combined, features = c("Kcng4"), cols = c("grey","blue"), order = F, min.cutoff = 0, reduction = "ref.umap", split.by = "condition2", keep.scale = "all")


#DimPlot(object =  multiome.combined, reduction = "ref.umap", group.by = "predicted.predicted_clusters", label=TRUE, pt.size = 0.25)



#renames clusters to add C in front, prevents a bug during DEG analysis
Idents(multiome.combined) <- "seurat_clusters"
clusterIDs <- paste0("C", multiome.combined$seurat_clusters) 
Idents(multiome.combined) <- clusterIDs
multiome.combined$renamedident <-Idents(multiome.combined)
DimPlot(object =  multiome.combined, reduction = "umap", label=TRUE, pt.size = 0.25)

# Differential gene expression 
#DefaultAssay(object = multiome.combined) <- "integrated"
Idents(multiome.combined) <-"condition"
#multiome.combined <- JoinLayers(object = multiome.combined)
OVAChase.response <- FindMarkers(multiome.combined, ident.1 = c("OVA_Chase"), ident.2 = c("PBS", "NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
#OVA.response2 <- FindMarkers(multiome.combined, ident.1 = c("OVA"), ident.2 = c("PBS", "NT"), slot = "data", test.use = "wilcox", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.1)

options(ggrepel.max.overlaps = Inf)
DefaultAssay(object = multiome.combined) <- "RNA"
aggregate_neurons <- AggregateExpression(multiome.combined, group.by = c("condition"), return.seurat = T, normalization.method = "LogNormalize", margin = 1, scale = 10000, slot = "data")
#aggregate_neurons <- AggregateExpression(multiome.combined, group.by = c("renamedident","condition"), return.seurat = T, normalization.method = "LogNormalize")
genes.to.label = c("Malat1")
#genes.to.label = rownames(head(OVA.response[order(-OVA.response$avg_log2FC),],n=10))

#top_genes <- OVA.response[OVA.response$avg_log2FC > 0, ]  # Filter for genes with ave_log2FC > 0
#top_genes <- top_genes[order(top_genes$p_val), ]  # Sort by p_val in ascending order
#genes.to.label = rownames(head(top_genes, n = 50))  # Select top 50 and return row names

# plotting
p1 <- CellScatter(aggregate_neurons, "OVA", "PBS", highlight = genes.to.label)
#p2 <- LabelPoints(plot = p1, points = genes.to.label,xnudge = -0.2, ynudge = 0.1, repel = T) 
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = T) 
p2 + xlim(0,6)+ylim(0,6) 

#p2 + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')

library(ggplot2)
library(ggrepel)
OVA.response$gene <- rownames(OVA.response)
ggplot(OVA.response, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() + ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val < 1e-8 & avg_log2FC >0.4, gene, "")), colour = "red", size = 3)
#ggplot(OVA.response2, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() + ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val < 1e-40 & avg_log2FC > 0, gene, "")), colour = "red", size = 3)

# find markers for clusters
new.markers <- FindAllMarkers(multiome.combined, assay = "RNA", features = NULL,  test.use = "roc")


# to map new data to reference dataset from 2020 Cell paper, aggregate.combined (AggregateData_190717_v2.RData)
aggregate.combined <- UpdateSeuratObject(aggregate.combined)
aggregate.combined <- RunUMAP(object = aggregate.combined, reduction = "pca", dims = 1:20, return.model=TRUE)
aggregate.combined$renamed_clusters <-Idents(aggregate.combined)
multiome_anchors <- FindTransferAnchors(reference = aggregate.combined, query = multiome.combined, dims = 1:30, reference.reduction = "pca")
multiome.combined <- MapQuery(anchorset = multiome_anchors, query = multiome.combined, reference = aggregate.combined, refdata = list(predicted_clusters = "renamed_clusters"), reference.reduction = "pca", reduction.model = "umap")

#DimPlot(multiome.combined, reduction = "ref.umap", label = T)
DimPlot(multiome.combined, group.by = "predicted.predicted_clusters",reduction = "ref.umap", label = T)

DimPlot(multiome.combined, group.by = "condition",reduction = "ref.umap", cols = c("blue", "orange", "red", "purple", "green"), pt.size = 0.5)

save(multiome.combined, file = "C:/Users/ekroj/MIT Dropbox/Emily Robinson/PrescottLab/Rcodes/multiome_combined_RNAseqAsthmaproject_250108.RData")

####
FeaturePlot(object= multiome.combined, features = c("Pvalb"), cols = c("grey","blue"), order = F, min.cutoff = 0, reduction = "ref.umap", split.by = "condition2", keep.scale = "all")


###

#determine which clusters have most gene expression changes in OVA conditions vs. PBS and NT
multiome.combined$cluster.stim <- paste(multiome.combined$predicted.predicted_clusters, multiome.combined$condition, sep = "_")
Idents(multiome.combined) <- multiome.combined$cluster.stim
JG1_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG1_OVA"), ident.2 = c("JG1_PBS", "JG1_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG2_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG2_OVA"), ident.2 = c("JG2_PBS", "JG2_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG3_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG3_OVA"), ident.2 = c("JG3_PBS", "JG3_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG4_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG4_OVA"), ident.2 = c("JG4_PBS", "JG4_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG5_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG5_OVA"), ident.2 = c("JG5_PBS", "JG5_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG6_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG6_OVA"), ident.2 = c("JG6_PBS", "JG6_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG7_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG7_OVA"), ident.2 = c("JG7_PBS", "JG7_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG8_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG8_OVA"), ident.2 = c("JG8_PBS", "JG8_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG9_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG9_OVA"), ident.2 = c("JG9_PBS", "JG9_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG10_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("JG10_OVA"), ident.2 = c("JG10_PBS", "JG10_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)

NG1_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG1_OVA"), ident.2 = c("NG1_PBS", "NG1_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG2_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG2_OVA"), ident.2 = c("NG2_PBS", "NG2_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG3_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG3_OVA"), ident.2 = c("NG3_PBS", "NG3_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG4_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG4_OVA"), ident.2 = c("NG4_PBS", "NG4_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG5_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG5_OVA"), ident.2 = c("NG5_PBS", "NG5_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG6_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG6_OVA"), ident.2 = c("NG6_PBS", "NG6_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG7_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG7_OVA"), ident.2 = c("NG7_PBS", "NG7_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG8_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG8_OVA"), ident.2 = c("NG8_PBS", "NG8_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG9_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG9_OVA"), ident.2 = c("NG9_PBS", "NG9_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG10_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG10_OVA"), ident.2 = c("NG10_PBS", "NG10_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG11_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG11_OVA"), ident.2 = c("NG11_PBS", "NG11_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG12_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG12_OVA"), ident.2 = c("NG12_PBS", "NG12_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG13_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG13_OVA"), ident.2 = c("NG13_PBS", "NG13_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
#NG14_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG14_OVA"), ident.2 = c("NG14_PBS", "NG14_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG15_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG15_OVA"), ident.2 = c("NG15_PBS", "NG15_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG16_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG16_OVA"), ident.2 = c("NG16_PBS", "NG16_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG17_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG17_OVA"), ident.2 = c("NG17_PBS", "NG17_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG18_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG18_OVA"), ident.2 = c("NG18_PBS", "NG18_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG19_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG19_OVA"), ident.2 = c("NG19_PBS", "NG19_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG20_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG20_OVA"), ident.2 = c("NG20_PBS", "NG20_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG21_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG21_OVA"), ident.2 = c("NG21_PBS", "NG21_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG22_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG22_OVA"), ident.2 = c("NG22_PBS", "NG22_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG23_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG23_OVA"), ident.2 = c("NG23_PBS", "NG23_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG24_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG24_OVA"), ident.2 = c("NG24_PBS", "NG24_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG25_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG25_OVA"), ident.2 = c("NG25_PBS", "NG25_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG26_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG26_OVA"), ident.2 = c("NG26_PBS", "NG26_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG27_OVA.degs <- FindMarkers(multiome.combined, ident.1 = c("NG27_OVA"), ident.2 = c("NG27_PBS", "NG27_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)

JG1_count <- nrow(JG1_OVA.degs[JG1_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG2_count <- nrow(JG2_OVA.degs[JG2_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG3_count <- nrow(JG3_OVA.degs[JG3_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG4_count <- nrow(JG4_OVA.degs[JG4_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG5_count <- nrow(JG5_OVA.degs[JG5_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG6_count <- nrow(JG6_OVA.degs[JG6_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG7_count <- nrow(JG7_OVA.degs[JG7_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG8_count <- nrow(JG8_OVA.degs[JG8_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG9_count <- nrow(JG9_OVA.degs[JG9_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
JG10_count <- nrow(JG10_OVA.degs[JG10_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0

NG1_count <- nrow(NG1_OVA.degs[NG1_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG2_count <- nrow(NG2_OVA.degs[NG2_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG3_count <- nrow(NG3_OVA.degs[NG3_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG4_count <- nrow(NG4_OVA.degs[NG4_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG5_count <- nrow(NG5_OVA.degs[NG5_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG6_count <- nrow(NG6_OVA.degs[NG6_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG7_count <- nrow(NG7_OVA.degs[NG7_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG8_count <- nrow(NG8_OVA.degs[NG8_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG9_count <- nrow(NG9_OVA.degs[NG9_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG10_count <- nrow(NG10_OVA.degs[NG10_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG11_count <- nrow(NG11_OVA.degs[NG11_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG12_count <- nrow(NG12_OVA.degs[NG12_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG13_count <- nrow(NG13_OVA.degs[NG13_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
#NG14_count <- nrow(NG14_OVA.degs[NG14_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG15_count <- nrow(NG15_OVA.degs[NG15_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG16_count <- nrow(NG16_OVA.degs[NG16_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG17_count <- nrow(NG17_OVA.degs[NG17_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG18_count <- nrow(NG18_OVA.degs[NG18_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG19_count <- nrow(NG19_OVA.degs[NG19_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG20_count <- nrow(NG20_OVA.degs[NG20_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG21_count <- nrow(NG21_OVA.degs[NG21_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG22_count <- nrow(NG22_OVA.degs[NG22_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG23_count <- nrow(NG23_OVA.degs[NG23_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG24_count <- nrow(NG24_OVA.degs[NG24_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG25_count <- nrow(NG25_OVA.degs[NG25_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG26_count <- nrow(NG26_OVA.degs[NG26_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0
NG27_count <- nrow(NG27_OVA.degs[NG27_OVA.degs$p_val_adj < 0.05,])  # Filter for genes with ave_log2FC > 0


#get count of differnetial genes from OVA Chase conditions
JG1_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG1_OVA_Chase"), ident.2 = c("JG1_PBS", "JG1_NT", "JG1_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG2_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG2_OVA_Chase"), ident.2 = c("JG2_PBS", "JG2_NT", "JG2_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG3_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG3_OVA_Chase"), ident.2 = c("JG3_PBS", "JG3_NT", "JG3_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG4_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG4_OVA_Chase"), ident.2 = c("JG4_PBS", "JG4_NT", "JG4_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG5_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG5_OVA_Chase"), ident.2 = c("JG5_PBS", "JG5_NT", "JG5_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
#JG6_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG6_OVA_Chase"), ident.2 = c("JG6_PBS", "JG6_NT", "JG6_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG6_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG6_OVA_Chase"), ident.2 = c("JG6_PBS", "JG6_NT"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG7_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG7_OVA_Chase"), ident.2 = c("JG7_PBS", "JG7_NT", "JG7_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG8_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG8_OVA_Chase"), ident.2 = c("JG8_PBS", "JG8_NT", "JG8_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG9_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG9_OVA_Chase"), ident.2 = c("JG9_PBS", "JG9_NT", "JG9_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
JG10_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("JG10_OVA_Chase"), ident.2 = c("JG10_PBS", "JG10_NT", "JG10_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)

NG1_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG1_OVA_Chase"), ident.2 = c("NG1_PBS", "NG1_NT", "NG1_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG2_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG2_OVA_Chase"), ident.2 = c("NG2_PBS", "NG2_NT", "NG2_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG3_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG3_OVA_Chase"), ident.2 = c("NG3_PBS", "NG3_NT", "NG3_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG4_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG4_OVA_Chase"), ident.2 = c("NG4_PBS", "NG4_NT", "NG4_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG5_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG5_OVA_Chase"), ident.2 = c("NG5_PBS", "NG5_NT", "NG5_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG6_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG6_OVA_Chase"), ident.2 = c("NG6_PBS", "NG6_NT", "NG6_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG7_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG7_OVA_Chase"), ident.2 = c("NG7_PBS", "NG7_NT", "NG7_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG8_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG8_OVA_Chase"), ident.2 = c("NG8_PBS", "NG8_NT", "NG8_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG9_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG9_OVA_Chase"), ident.2 = c("NG9_PBS", "NG9_NT", "NG9_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG10_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG10_OVA_Chase"), ident.2 = c("NG10_PBS", "NG10_NT", "NG10_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG11_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG11_OVA_Chase"), ident.2 = c("NG11_PBS", "NG11_NT", "NG11_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG12_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG12_OVA_Chase"), ident.2 = c("NG12_PBS", "NG12_NT", "NG12_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG13_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG13_OVA_Chase"), ident.2 = c("NG13_PBS", "NG13_NT", "NG13_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
#NG14_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG14_OVA_Chase"), ident.2 = c("NG14_PBS", "NG14_NT", "NG14_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG15_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG15_OVA_Chase"), ident.2 = c("NG15_PBS", "NG15_NT", "NG15_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG16_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG16_OVA_Chase"), ident.2 = c("NG16_PBS", "NG16_NT", "NG16_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG17_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG17_OVA_Chase"), ident.2 = c("NG17_PBS", "NG17_NT", "NG17_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG18_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG18_OVA_Chase"), ident.2 = c("NG18_PBS", "NG18_NT", "NG18_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG19_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG19_OVA_Chase"), ident.2 = c("NG19_PBS", "NG19_NT", "NG19_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG20_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG20_OVA_Chase"), ident.2 = c("NG20_PBS", "NG20_NT", "NG20_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG21_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG21_OVA_Chase"), ident.2 = c("NG21_PBS", "NG21_NT", "NG21_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG22_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG22_OVA_Chase"), ident.2 = c("NG22_PBS", "NG22_NT", "NG22_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG23_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG23_OVA_Chase"), ident.2 = c("NG23_PBS", "NG23_NT", "NG23_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG24_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG24_OVA_Chase"), ident.2 = c("NG24_PBS", "NG24_NT", "NG24_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG25_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG25_OVA_Chase"), ident.2 = c("NG25_PBS", "NG25_NT", "NG25_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG26_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG26_OVA_Chase"), ident.2 = c("NG26_PBS", "NG26_NT", "NG26_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)
NG27_OVAChase.degs <- FindMarkers(multiome.combined, ident.1 = c("NG27_OVA_Chase"), ident.2 = c("NG27_PBS", "NG27_NT", "NG27_PBS_Chase"), slot = "counts", test.use = "DESeq2", verbose = FALSE, only.pos = F, logfc.threshold = 0.5,  min.pct = 0.5)


JG1_Chase_count <- nrow(JG1_OVAChase.degs[JG1_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG2_Chase_count <- nrow(JG1_OVAChase.degs[JG2_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG3_Chase_count <- nrow(JG5_OVAChase.degs[JG3_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG4_Chase_count <- nrow(JG5_OVAChase.degs[JG4_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG5_Chase_count <- nrow(JG5_OVAChase.degs[JG2_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG6_Chase_count <- nrow(JG6_OVAChase.degs[JG6_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG7_Chase_count <- nrow(JG7_OVAChase.degs[JG7_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG8_Chase_count <- nrow(JG8_OVAChase.degs[JG8_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG9_Chase_count <- nrow(JG9_OVAChase.degs[JG9_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
JG10_Chase_count <- nrow(JG10_OVAChase.degs[JG10_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0

NG1_Chase_count <- nrow(NG1_OVAChase.degs[NG1_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG2_Chase_count <- nrow(NG2_OVAChase.degs[NG2_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG3_Chase_count <- nrow(NG3_OVAChase.degs[NG3_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG4_Chase_count <- nrow(NG5_OVAChase.degs[NG4_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG5_Chase_count <- nrow(NG5_OVAChase.degs[NG2_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG6_Chase_count <- nrow(NG6_OVAChase.degs[NG6_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG7_Chase_count <- nrow(NG7_OVAChase.degs[NG7_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG8_Chase_count <- nrow(NG8_OVAChase.degs[NG8_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG9_Chase_count <- nrow(NG9_OVAChase.degs[NG9_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG10_Chase_count <- nrow(NG10_OVAChase.degs[NG10_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG11_Chase_count <- nrow(NG11_OVAChase.degs[NG11_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG12_Chase_count <- nrow(NG12_OVAChase.degs[NG12_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG13_Chase_count <- nrow(NG13_OVAChase.degs[NG13_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
#NG14_Chase_count <- nrow(NG15_OVAChase.degs[NG14_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG15_Chase_count <- nrow(NG15_OVAChase.degs[NG15_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG16_Chase_count <- nrow(NG16_OVAChase.degs[NG16_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG17_Chase_count <- nrow(NG17_OVAChase.degs[NG17_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG18_Chase_count <- nrow(NG18_OVAChase.degs[NG18_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG19_Chase_count <- nrow(NG19_OVAChase.degs[NG19_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG20_Chase_count <- nrow(NG20_OVAChase.degs[NG20_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG21_Chase_count <- nrow(NG21_OVAChase.degs[NG21_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG22_Chase_count <- nrow(NG22_OVAChase.degs[NG22_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG23_Chase_count <- nrow(NG23_OVAChase.degs[NG23_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG24_Chase_count <- nrow(NG24_OVAChase.degs[NG24_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG25_Chase_count <- nrow(NG25_OVAChase.degs[NG25_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG26_Chase_count <- nrow(NG26_OVAChase.degs[NG26_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
NG27_Chase_count <- nrow(NG27_OVAChase.degs[NG27_OVAChase.degs$p_val_adj < 0.05,])  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0


#write.csv(NG11_OVA.degs,"~/Downloads/NG11_OVA.degs.csv", row.names = T)
write.csv(NG11_OVA.degs[NG11_OVA.degs$p_val_adj < 0.05,],"C:/Users/ekroj/OneDrive/Desktop/RPlots_241031/NG11_OVA.degs.csv", row.names = T)
write.csv(NG11_OVAChase.degs[NG11_OVAChase.degs$p_val_adj < 0.05,],"C:/Users/ekroj/OneDrive/Desktop/RPlots_241031/NG11_OVAChase.degs.csv", row.names = T)


#make a venn diagram of gene overlaps
install.packages("ggVennDiagram")
library(ggVennDiagram)
NG11_Venn <- list(A = rownames(NG11_OVA.degs[NG11_OVA.degs$p_val_adj < 0.05,]), B = rownames(NG11_OVAChase.degs[NG11_OVAChase.degs$p_val_adj < 0.05,]))
NG13_Venn <- list(A = rownames(NG13_OVA.degs[NG13_OVA.degs$p_val_adj < 0.05,]), B = rownames(NG13_OVAChase.degs[NG13_OVAChase.degs$p_val_adj < 0.05,]))
ggVennDiagram(NG11_Venn, category.names = c("OVA","OVA Chase"),label_alpha = 0)+coord_flip() + theme(legend.position = "none") + scale_fill_gradient(low = "#CC3300", high = "#FFFCCC")
ggVennDiagram(NG13_Venn, category.names = c("OVA","OVA Chase"),label_alpha = 0)+coord_flip() + theme(legend.position = "none") + scale_fill_gradient(low = "#CC3300", high = "#FFFCCC")


#To find list of shared genes between OVA and OVA Chase
ng11_genes <- rownames(NG11_OVA.degs[NG11_OVA.degs$p_val_adj < 0.05, ])
ng11_chase_genes <- rownames(NG11_OVAChase.degs[NG11_OVAChase.degs$p_val_adj < 0.05, ])
shared_genes_ng11 <- intersect(ng11_genes, ng11_chase_genes)
print(shared_genes_ng11)     


######### My work ########

# Convert row names to a column for NG11_OVA.degs
NG11_OVA.degs <- as.data.frame(NG11_OVA.degs)  # Convert to data frame if not already
NG11_OVA.degs$gene <- rownames(NG11_OVA.degs)

# Convert row names to a column for NG11_OVAChase.degs
NG11_OVAChase.degs <- as.data.frame(NG11_OVAChase.degs)  # Convert to data frame if not already
NG11_OVAChase.degs$gene <- rownames(NG11_OVAChase.degs)



# Step 2: Identify Shared Significant Genes
shared_genes <- merge(
  NG11_OVA.degs[NG11_OVA.degs$p_val_adj < 0.05, c("avg_log2FC", "p_val_adj")],
  NG11_OVAChase.degs[NG11_OVAChase.degs$p_val_adj < 0.05, c("avg_log2FC", "p_val_adj")],
  by = "row.names",  # Merge by gene name (row name)
  suffixes = c("_OVA", "_OVAChase")
)

# Convert to data frame and set the column names correctly
shared_genes_df <- as.data.frame(shared_genes)
colnames(shared_genes_df) <- c("gene", "avg_log2FC_OVA", "p_val_adj_OVA", "avg_log2FC_OVAChase", "p_val_adj_OVAChase")

# Step 3: Filter for Positive Log2FC Genes
#positive_shared_genes <- shared_genes_df[
#  shared_genes_df$avg_log2FC_OVA > 0 & shared_genes_df$avg_log2FC_OVAChase > 0, 
#]

# Step 4: Prepare Data for Dot Plot
dot_plot_data <- data.frame(
  gene = shared_genes_df$gene,
  log2FC_OVA = shared_genes_df$avg_log2FC_OVA,
  log2FC_OVAChase = shared_genes_df$avg_log2FC_OVAChase,
  p_val_OVA = -log10(shared_genes_df$p_val_adj_OVA),
  p_val_OVAChase = -log10(shared_genes_df$p_val_adj_OVAChase)
)

# Step 4: Create the Dot Plot with Arrows and Gene Names
library(ggplot2)

# Plotting
ggplot(dot_plot_data) +
  # Arrows connecting the points
  geom_segment(aes(x = p_val_OVA, y = log2FC_OVA, 
                   xend = p_val_OVAChase, yend = log2FC_OVAChase), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "gray", size = 0.5) +
  # Points for OVA
  geom_point(aes(x = p_val_OVA, y = log2FC_OVA, color = "OVA"), size = 3) +
  # Points for OVA Chase
  geom_point(aes(x = p_val_OVAChase, y = log2FC_OVAChase, color = "OVA Chase"), size = 3) +
  # Gene names along the arrows
  geom_text(aes(x = (p_val_OVA + p_val_OVAChase) / 2, 
                y = (log2FC_OVA + log2FC_OVAChase) / 2, 
                label = gene), 
            size = 3, vjust = -0.5, color = "black") +
  labs(
    title = "Dot Plot of Shared Significant Genes",
    subtitle = paste("Total shared significant genes:", nrow(dot_plot_data)),
    x = "-log10(p-value)",
    y = "Log2 Fold Change"
  ) +
  scale_color_manual(values = c("OVA" = "blue", "OVA Chase" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank()) +  # Remove legend title
  guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust legend point size



# Step 1: Calculate Absolute Log2 Fold Change Differences
shared_genes_df$abs_log2FC_OVA <- abs(shared_genes_df$avg_log2FC_OVA)
shared_genes_df$abs_log2FC_OVAChase <- abs(shared_genes_df$avg_log2FC_OVAChase)

# Step 2: Calculate the difference in absolute values
shared_genes_df$abs_diff <- abs(shared_genes_df$abs_log2FC_OVA - shared_genes_df$abs_log2FC_OVAChase)

# Step 3: Filter for genes that:
# - Increase in absolute log2FC value in OVAChase
# - Have abs log2FC > 0.5 in the OVA case
filtered_genes <- shared_genes_df[
  (shared_genes_df$abs_log2FC_OVAChase > shared_genes_df$abs_log2FC_OVA) & 
    (shared_genes_df$abs_log2FC_OVA > 0.5), 
]

# Now filtered_genes contains the genes that meet both criteria

# Prepare Data for Dot Plot
dot_plot_data <- data.frame(
  gene = filtered_genes$gene,
  log2FC_OVA = filtered_genes$avg_log2FC_OVA,
  log2FC_OVAChase = filtered_genes$avg_log2FC_OVAChase,
  p_val_OVA = -log10(filtered_genes$p_val_adj_OVA),
  p_val_OVAChase = -log10(filtered_genes$p_val_adj_OVAChase)
)

# Continue with plotting as before


# Plotting
ggplot(dot_plot_data) +
  # Arrows connecting the points
  geom_segment(aes(x = p_val_OVA, y = log2FC_OVA, 
                   xend = p_val_OVAChase, yend = log2FC_OVAChase), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "gray", size = 0.5) +
  # Points for OVA
  geom_point(aes(x = p_val_OVA, y = log2FC_OVA, color = "OVA"), size = 3) +
  # Points for OVA Chase
  geom_point(aes(x = p_val_OVAChase, y = log2FC_OVAChase, color = "OVA Chase"), size = 3) +
  # Gene names along the arrows
  geom_text(aes(x = (p_val_OVA + p_val_OVAChase) / 2, 
                y = (log2FC_OVA + log2FC_OVAChase) / 2, 
                label = gene), 
            size = 3, vjust = -0.5, color = "black") +
  labs(
    title = "Dot Plot of Shared Significant Genes",
    subtitle = paste("Total shared significant genes:", nrow(dot_plot_data)),
    x = "-log10(p-value)",
    y = "Log2 Fold Change"
  ) +
  scale_color_manual(values = c("OVA" = "blue", "OVA Chase" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank()) +  # Remove legend title
  guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust legend point size






##Making a dot plot with interesting genes 


NG11_NT_markers <- FindMarkers(multiome.combined, 
                               ident.1 = c("NG11_NT", "NG11_PBS"),  # Directly use the cluster name
                               ident.2 = c("NG1_NT", "NG2_NT", "NG3_NT", "NG4_NT", "NG5_NT", 
                                           "NG6_NT", "NG7_NT", "NG8_NT", "NG9_NT", "NG10_NT", 
                                           "NG12_NT", "NG13_NT", "NG14_NT", "NG15_NT", "NG16_NT", 
                                           "NG17_NT", "NG18_NT", "NG19_NT", "NG20_NT", "NG21_NT", 
                                           "NG22_NT", "NG23_NT", "NG24_NT", "NG25_NT", "NG26_NT", 
                                           "NG27_NT", "JG1_NT", "JG2_NT", "JG3_NT", "JG4_NT", 
                                           "JG5_NT", "JG6_NT", "JG7_NT", "JG8_NT", "JG9_NT", 
                                           "JG10_NT", "NG1_PBS", "NG2_PBS", "NG3_PBS", "NG4_PBS", 
                                           "NG5_PBS", "NG6_PBS", "NG7_PBS", "NG8_PBS", "NG9_PBS", 
                                           "NG10_PBS", "NG12_PBS", "NG13_PBS", 
                                           "NG15_PBS", "NG16_PBS", "NG17_PBS", "NG18_PBS", 
                                           "NG19_PBS", "NG20_PBS", "NG21_PBS", "NG22_PBS", 
                                           "NG23_PBS", "NG24_PBS", "NG25_PBS", "NG26_PBS", 
                                           "NG27_PBS", "JG1_PBS", "JG2_PBS", "JG3_PBS", 
                                           "JG4_PBS", "JG5_PBS", "JG6_PBS", "JG7_PBS", 
                                           "JG8_PBS", "JG9_PBS", "JG10_PBS"),  # List all other clusters directly
                               slot = "counts", 
                               test.use = "DESeq2", 
                               verbose = FALSE)


head(NG11_NT_count)

NG11_NT_count <- NG11_NT_markers[NG11_NT_markers$p_val_adj < 0.05,]
NG11_NT_count <- NG11_NT_count[NG11_NT_markers$avg_log2FC > 0.2,]  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0
#NG11_NT_count <- NG11_NT_count[NG11_NT_markers$pct.2 < 0.4,]  # Filter for _NT", "X_PBSChase" with ave_log2FC > 0


# Step 1: Subset the Seurat object to include only NT and PBS conditions
multiome_subset <- subset(multiome.combined, condition %in% c("NT", "PBS"))


# Step 2: Ensure the correct column exists for subtypes. Let's assume it's 'subtype' or another name
# Check if a column that defines subtypes already exists:
if (!"subtype" %in% colnames(multiome_subset@meta.data)) {
  # If not, you need to create this column (this step depends on how you define subtypes)
  multiome_subset$subtype <- as.character(multiome_subset$predicted.predicted_clusters)  # Change logic as needed
}

# Step 3: Create the 'NG11_vs_rest_subtypes' column to include all unique subtypes
multiome_subset$NG11_vs_rest_subtypes <- factor(multiome_subset$subtype)

# Step 4: Check the unique values in the 'NG11_vs_rest_subtypes' column
unique(multiome_subset$NG11_vs_rest_subtypes)

# This should show all your subtypes individually without collapsing NG11 against Rest.

# Step 3: Define the genes you're interested in
genes <- c("Cysltr2", "Dbh", "C1ql2", "Piezo2", "Npy1r", "P2ry1", "Gabra1", "F2r", "S1pr3", 
           "Crhr2", "Gpr149", "Gpr174", "Gpr176", "Npy2r", "Drd2", "Ptger4", 
           "Trpa1", "Trpv1", "Tac1", "Nos1", "Calca", "Phox2b", "Prdm12")

# Extract the normalized gene expression data from the 'data' slot of the 'RNA' assay
# Step 1: Access the 'counts' matrix (raw data)
#rna_data <- as.matrix(multiome_subset[["RNA"]]@layers$scale.data)


# Step 2: Clip values to the range of 0 to 2

#rna_data[rna_data < 0] <- 0

# Step 3: Assign the modified data back to the Seurat object
#multiome_subset[["RNA"]]@layers$scale.data <- rna_data

# Step 4: Generate the DotPlot with the clipped raw counts
# Check the values before plotting

# Test a custom color scale (optional)


DotPlot(object = multiome_subset, 
        features = genes, 
        assay = "RNA", 
        cols = c("gold", "blue"), 
        scale.by = "size", 
        group.by = "NG11_vs_rest_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ng11_genes <- rownames(NG11_OVA.degs[NG11_OVA.degs$p_val_adj < 0.05, ])
ng11_chase_genes <- rownames(NG11_OVAChase.degs[NG11_OVAChase.degs$p_val_adj < 0.05, ])
shared_genes_ng11_all <- intersect(ng11_genes, ng11_chase_genes, NG11_NT_count)


print(shared_genes_ng11_all)     



head(NG11_NT_count)
write.csv(NG11_NT_count,"C:/Users/ekroj/OneDrive/Desktop/RPlots_241031/NG11_NT_counts.csv", row.names = T)



library(ggplot2)

# Define the genes for each condition
control_genes <- c("Gene1", "Gene2", "Gene3")  # Replace with actual control genes
ova_genes <- c("Npy1r", "Asic3", "Il4ra")  # Example OVA genes
ova_chase_genes <- c("Tlr4", "Ptger1", "Il13ra1")  # Example OVA Chase genes

# Combine all genes into a single vector for plotting
all_genes <- c(control_genes, ova_genes, ova_chase_genes)


# Subset for NG11 and all other clusters
ng11_vs_others <- subset(multiome.combined, 
                         subset = Idents(new.combined) %in% c("NG11", "NG1", "NG2", "NG3", "NG4", "NG5", 
                                                              "NG6", "NG7", "NG8", "NG9", "NG10", 
                                                              "NG12", "NG13", "NG14", "NG15", "NG16", 
                                                              "NG17", "NG18", "NG19", "NG20", "NG21", 
                                                              "NG22", "NG23", "NG24", "NG25", "NG26", 
                                                              "NG27", "JG1", "JG2", "JG3", "JG4", "JG5", 
                                                              "JG6", "JG7", "JG8", "JG9", "JG10"))

ng11_vs_others@active.ident <- factor(ng11_vs_others@active.ident, 
                                      levels = c("NG11", "NG1", "NG2", "NG3", "NG4", "NG5", 
                                                 "NG6", "NG7", "NG8", "NG9", "NG10", 
                                                 "NG12", "NG13", "NG14", "NG15", "NG16", 
                                                 "NG17", "NG18", "NG19", "NG20", "NG21", 
                                                 "NG22", "NG23", "NG24", "NG25", "NG26", 
                                                 "NG27", "JG1", "JG2", "JG3", "JG4", "JG5", 
                                                 "JG6", "JG7", "JG8", "JG9", "JG10"))


# Create the Dot Plot
dot_plot <- DotPlot(object = ng11_vs_others, features = all_genes, assay = "RNA", 
                    col.min = 0, col.max = 3, scale.by = 'size') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Dot Plot of NG11 vs All Other Subtypes")



# Add facetting to separate conditions visually (if applicable)
dot_plot <- dot_plot + 
  facet_wrap(~condition, scales = "free")  # Ensure 'condition' variable is in your data

# Display the plot
print(dot_plot)


# Make a volcano plot for NG11, NG13 and JG1

# Load necessary libraries
library(ggplot2)
library(ggrepel)  # Optional: for better label placement

# Assuming NG11_OVA.degs is your dataframe
# Calculate -log10(p-value) for better visualization
NG11_OVA.degs$neg_log10_pval <- -log10(NG11_OVA.degs$p_val_adj)

# Create a new column for significant genes based on the criteria
# Adjust these thresholds (p-value < 0.05 and |log2FC| > 1) based on your specific needs
NG11_OVA.degs$significant <- ifelse(NG11_OVA.degs$p_val_adj < 0.05 & abs(NG11_OVA.degs$avg_log2FC) > 0.5, "Significant", "Not Significant")

# Subset significant genes
significant_genes <- subset(NG11_OVA.degs, significant == "Significant")

# Volcano plot with labels for significant genes
ggplot(NG11_OVA.degs, aes(x = avg_log2FC, y = neg_log10_pval, color = significant)) +
  geom_point(size = 2) +  # Points for all genes
  scale_color_manual(values = c("gray", "red")) +  # Red for significant genes, gray for others
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Volcano Plot for Upregulated Genes in NG11",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted p-value)"
  ) +
  # Add labels for significant genes using rownames() on the subset
  geom_text_repel(
    data = significant_genes,
    aes(label = rownames(significant_genes)),  # Use rownames() from the subset
    size = 3,                              # Adjust text size
    box.padding = 0.5,                      # Adjust padding around labels
    point.padding = 0.5,                    # Adjust the space between points and labels
    color = "black",                        # Label color
    max.overlaps = 100 
  )


# NG11 OVA Chase


# Load necessary libraries
library(ggplot2)
library(ggrepel)  # Optional: for better label placement

# Assuming NG11_OVA.degs is your dataframe
# Calculate -log10(p-value) for better visualization
NG11_OVAChase.degs$neg_log10_pval <- -log10(NG11_OVAChase.degs$p_val_adj)

# Create a new column for significant genes based on the criteria
# Adjust these thresholds (p-value < 0.05 and |log2FC| > 1) based on your specific needs
NG11_OVAChase.degs$significant <- ifelse(NG11_OVAChase.degs$p_val_adj < 0.05 & abs(NG11_OVAChase.degs$avg_log2FC) > 0.5, "Significant", "Not Significant")

# Subset significant genes
significant_genes <- subset(NG11_OVAChase.degs, significant == "Significant")

# Volcano plot with labels for significant genes
ggplot(NG11_OVAChase.degs, aes(x = avg_log2FC, y = neg_log10_pval, color = significant)) +
  geom_point(size = 2) +  # Points for all genes
  scale_color_manual(values = c("gray", "red")) +  # Red for significant genes, gray for others
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Volcano Plot for Upregulated Genes in NG11 in OVA Chase",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted p-value)"
  ) +
  # Add labels for significant genes using rownames() on the subset
  geom_text_repel(
    data = significant_genes,
    aes(label = rownames(significant_genes)),  # Use rownames() from the subset
    size = 3,                              # Adjust text size
    box.padding = 0.5,                      # Adjust padding around labels
    point.padding = 0.5,                    # Adjust the space between points and labels
    color = "black",                        # Label color
    max.overlaps = 100 
  )

# NG13 OVA Volcano Plot

# Load necessary libraries
library(ggplot2)
library(ggrepel)  # Optional: for better label placement

# Assuming NG11_OVA.degs is your dataframe
# Calculate -log10(p-value) for better visualization
NG13_OVA.degs$neg_log10_pval <- -log10(NG13_OVA.degs$p_val_adj)

# Create a new column for significant genes based on the criteria
# Adjust these thresholds (p-value < 0.05 and |log2FC| > 1) based on your specific needs
NG13_OVA.degs$significant <- ifelse(NG13_OVA.degs$p_val_adj < 0.05 & abs(NG13_OVA.degs$avg_log2FC) > 0.5, "Significant", "Not Significant")

# Subset significant genes
significant_genes <- subset(NG13_OVA.degs, significant == "Significant")

# Volcano plot with labels for significant genes
ggplot(NG13_OVA.degs, aes(x = avg_log2FC, y = neg_log10_pval, color = significant)) +
  geom_point(size = 2) +  # Points for all genes
  scale_color_manual(values = c("gray", "red")) +  # Red for significant genes, gray for others
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Volcano Plot for Upregulated Genes in NG13 in OVA Condition",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted p-value)"
  ) +
  # Add labels for significant genes using rownames() on the subset
  geom_text_repel(
    data = significant_genes,
    aes(label = rownames(significant_genes)),  # Use rownames() from the subset
    size = 3,                              # Adjust text size
    box.padding = 0.5,                      # Adjust padding around labels
    point.padding = 0.5,                    # Adjust the space between points and labels
    color = "black",                        # Label color
    max.overlaps = 100 
  )


# JG1 OVA Volcano Plot

# Load necessary libraries
library(ggplot2)
library(ggrepel)  # Optional: for better label placement

# Assuming NG11_OVA.degs is your dataframe
# Calculate -log10(p-value) for better visualization
JG1_OVA.degs$neg_log10_pval <- -log10(JG1_OVA.degs$p_val_adj)

# Create a new column for significant genes based on the criteria
# Adjust these thresholds (p-value < 0.05 and |log2FC| > 1) based on your specific needs
JG1_OVA.degs$significant <- ifelse(JG1_OVA.degs$p_val_adj < 0.05 & abs(JG1_OVA.degs$avg_log2FC) > 0.5, "Significant", "Not Significant")

# Subset significant genes
significant_genes <- subset(JG1_OVA.degs, significant == "Significant")

# Volcano plot with labels for significant genes
ggplot(JG1_OVA.degs, aes(x = avg_log2FC, y = neg_log10_pval, color = significant)) +
  geom_point(size = 2) +  # Points for all genes
  scale_color_manual(values = c("gray", "red")) +  # Red for significant genes, gray for others
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Volcano Plot for Upregulated Genes in JG1 in OVA Condition",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted p-value)"
  ) +
  # Add labels for significant genes using rownames() on the subset
  geom_text_repel(
    data = significant_genes,
    aes(label = rownames(significant_genes)),  # Use rownames() from the subset
    size = 3,                              # Adjust text size
    box.padding = 0.5,                      # Adjust padding around labels
    point.padding = 0.5,                    # Adjust the space between points and labels
    color = "black",                        # Label color
    max.overlaps = 100 
  )



# Making Violin Plot for GPCRs in NG11


# Load necessary libraries
library(ggplot2)
library(ggsignif)

# Subset the Seurat object to include only the desired conditions
multiome_subset <- subset(multiome_subset, condition %in% c("NT", "PBS", "OVA", "OVA_Chase"))

# Reorder the factor levels of the 'condition' column in your metadata
multiome_subset$condition <- factor(multiome_subset$condition, 
                                    levels = c("NT", "PBS", "OVA", "OVA_Chase"))

# Fetch expression data for "gene of interest"
expr_data <- FetchData(multiome_subset, vars = "F2r")

# Extract the condition information from the Seurat object
condition_data <- multiome_subset$condition

# Create a data frame with both condition and expression data
expr_data_with_condition <- data.frame(
  F2r = expr_data$F2r,
  condition = condition_data
)

# Perform pairwise t-tests between conditions
t_test_PBS_OVA <- t.test(expr_data_with_condition$F2r[expr_data_with_condition$condition == "PBS"], 
                         expr_data_with_condition$F2r[expr_data_with_condition$condition == "OVA"])

t_test_PBS_OVA_Chase <- t.test(expr_data_with_condition$F2r[expr_data_with_condition$condition == "PBS"], 
                               expr_data_with_condition$F2r[expr_data_with_condition$condition == "OVA_Chase"])

t_test_OVA_OVA_Chase <- t.test(expr_data_with_condition$F2r[expr_data_with_condition$condition == "OVA"], 
                               expr_data_with_condition$F2r[expr_data_with_condition$condition == "OVA_Chase"])

# Function to convert p-value to significance stars
get_significance <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")  # Not significant
  }
}

# Create a vector of significance annotations based on p-values
annotations <- c(
  get_significance(t_test_PBS_OVA$p.value),
  get_significance(t_test_PBS_OVA_Chase$p.value),
  get_significance(t_test_OVA_OVA_Chase$p.value)
)

# Create the violin plot for the gene of interest
p <- ggplot(expr_data_with_condition, aes(x = condition, y = F2r, fill = condition)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Define y_positions for each comparison to avoid overlap
y_positions <- c(5, 6, 7)  # Adjust these values depending on the scale of your plot

# Add significance bars with annotations and custom y_position
p + geom_signif(
  comparisons = list(c("PBS", "OVA"), c("PBS", "OVA_Chase"), c("OVA", "OVA_Chase")),  # Specify condition pairs to compare
  annotations = annotations,  # Add annotations as significance stars
  map_signif_level = TRUE,     # Show p-value as stars
  test = "t.test",             # Perform t-test automatically
  y_position = y_positions     # Specify y_position for each comparison
) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
  axis.text.y = element_text(size = 14),  # Increase y-axis text size
  axis.title.x = element_text(size = 16),  # Increase x-axis title size
  axis.title.y = element_text(size = 16),  # Increase y-axis title size
  plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size
  legend.text = element_text(size = 14),  # Increase legend text size
  annotation.text = element_text(size = 14)  # Increase annotation (significance) text size
)



# Dot plot for all subtypes


# define genes of interest

#"Cysltr2", "Dbh", "C1ql2", "Piezo2", "Npy1r", "P2ry1", "Gabra1",

genes <- c("F2r", "S1pr3", 
           "Crhr2", "Gpr149", "Gpr174", "Npy2r", "Drd2", "Ptger4", "Il4ra", "Tlr4",
           "Trpa1", "Trpv1", "Tac1", "Nos1", "Vip", "Calca", "Phox2b", "Prdm12")

# Create the Dot Plot
DotPlot(object = multiome.combined, features = genes, assay = "RNA", 
        +                     col.min = 0, col.max = 3, scale = T, scale.by = 'size', dot.scale = 4,  group.by = "predicted.predicted_clusters") + 
  +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  +   labs(title = "Dot Plot of All Subtypes")


