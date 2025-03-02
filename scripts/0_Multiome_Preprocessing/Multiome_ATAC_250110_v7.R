library(Seurat)
library(Signac)
library(GenomicFeatures)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(Rsamtools)


# NT: the 10x file contains both data types. 
NT.data <- Read10X(data.dir = "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/1_NT_filtered_feature_bc_matrix")

# extract RNA and ATAC data
rna_counts <- NT.data$`Gene Expression`
atac_counts <- NT.data$Peaks

# Create Seurat object
NT <- CreateSeuratObject(counts = rna_counts)
NT[["percent.mt"]] <- PercentageFeatureSet(NT, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/25_01_05_ATAC/NT/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
NT[["ATAC"]] <- chrom_assay


# PBS
# the 10x file contains both data types. 
PBS.data <- Read10X(data.dir = "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/2_PBS_filtered_feature_bc_matrix")
# extract RNA and ATAC data
rna_counts <- PBS.data$`Gene Expression`
atac_counts <- PBS.data$Peaks
# Create Seurat object
PBS <- CreateSeuratObject(counts = rna_counts)
PBS[["percent.mt"]] <- PercentageFeatureSet(PBS, pattern = "^MT-")
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
frag.file <- "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/25_01_05_ATAC/PBS/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
PBS[["ATAC"]] <- chrom_assay


# OVA
# the 10x file contains both data types. 
OVA.data <- Read10X(data.dir = "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/3_OVA_filtered_feature_bc_matrix")
# extract RNA and ATAC data
rna_counts <- OVA.data$`Gene Expression`
atac_counts <- OVA.data$Peaks
# Create Seurat object
OVA <- CreateSeuratObject(counts = rna_counts)
OVA[["percent.mt"]] <- PercentageFeatureSet(OVA, pattern = "^MT-")
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
frag.file <- "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/25_01_05_ATAC/OVA/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
OVA[["ATAC"]] <- chrom_assay


# the 10x file contains both data types. 
PBS_C.data <- Read10X(data.dir = "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/4_PBS_C_filtered_feature_bc_matrix")
# extract RNA and ATAC data
rna_counts <- PBS_C.data$`Gene Expression`
atac_counts <- PBS_C.data$Peaks
# Create Seurat object
PBS_C <- CreateSeuratObject(counts = rna_counts)
PBS_C[["percent.mt"]] <- PercentageFeatureSet(PBS_C, pattern = "^MT-")
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
frag.file <- "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/25_01_05_ATAC/PBS_C/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
PBS_C[["ATAC"]] <- chrom_assay


# OVA Chase
# the 10x file contains both data types. 
OVA_C.data <- Read10X(data.dir = "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/24_09_27/5_OVA_C_filtered_feature_bc_matrix")
# extract RNA and ATAC data
rna_counts <- OVA_C.data$`Gene Expression`
atac_counts <- OVA_C.data$Peaks
# Create Seurat object
OVA_C <- CreateSeuratObject(counts = rna_counts)
OVA_C[["percent.mt"]] <- PercentageFeatureSet(OVA_C, pattern = "^MT-")
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
frag.file <- "C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Results/scRNAseqanalysis/25_01_05_ATAC/OVA_C/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
OVA_C[["ATAC"]] <- chrom_assay


#Merge all seurat objects - this step requires A LOT of RAM, don't do too much on your computer while running this, can cause crash
multiome.combined <- merge(NT, y = list(PBS, OVA, PBS_C, OVA_C), add.cell.ids = c("NT", "PBS", "OVA", "PBS_C", "OVA_C"))

#save(multiome.combined, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Rcodes/CombinedMultiomeData_250110.RData")



#Subsetting out "real cells"
VlnPlot(multiome.combined, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
multiome.combined <- subset(
  x = multiome.combined,
  subset = nCount_ATAC >= 1000)

# RNA analysis
DefaultAssay(multiome.combined) <- "RNA"
multiome.combined <- NormalizeData(multiome.combined)
multiome.combined <- FindVariableFeatures(multiome.combined)
multiome.combined <- ScaleData(multiome.combined)
multiome.combined <- RunPCA(multiome.combined)
multiome.combined <- RunUMAP(multiome.combined, dims = 1:30)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(multiome.combined) <- "ATAC"
multiome.combined  <- RunTFIDF(multiome.combined )
multiome.combined  <- FindTopFeatures(multiome.combined , min.cutoff = 'q0')
multiome.combined  <- RunSVD(multiome.combined)
multiome.combined  <- RunUMAP(multiome.combined , reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(multiome.combined, group.by = "orig.ident", label = FALSE)

# quantify gene activity
DefaultAssay(multiome.combined) <- "ATAC"
gene.activities <- GeneActivity(multiome.combined)

# add gene activities as a new assay
multiome.combined[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

#Subset out only cells that are neurons
DefaultAssay(multiome.combined) <- "RNA"
FeaturePlot(multiome.combined, features = "Slc17a6", reduction = "umap.atac", order = TRUE)
multiome.combined$is_neuron <- ifelse(
  LayerData(multiome.combined[["RNA"]], layer = "counts.1")["Slc17a6",] > 0,
  "Neuron",
  "Non-neuron"
)

multiome.combined$is_neuron <- ifelse(
  LayerData(multiome.combined[["RNA"]], layer = "counts.2")["Slc17a6",] > 0,
  "Neuron",
  "Non-neuron"
)


multiome.combined$is_neuron <- ifelse(
  LayerData(multiome.combined[["RNA"]], layer = "counts.3")["Slc17a6",] > 0,
  "Neuron",
  "Non-neuron"
)

multiome.combined$is_neuron <- ifelse(
  LayerData(multiome.combined[["RNA"]], layer = "counts.4")["Slc17a6",] > 0,
  "Neuron",
  "Non-neuron"
)
multiome.combined$is_neuron <- ifelse(
  LayerData(multiome.combined[["RNA"]], layer = "counts.5")["Slc17a6",] > 0,
  "Neuron",
  "Non-neuron"
)

DimPlot(multiome.combined, reduction = "umap.atac", group.by = "is_neuron")
multiome.combined.neurons <- subset(multiome.combined, subset = is_neuron == "Neuron")
DimPlot(multiome.combined.neurons, reduction = "umap.atac")
multiome.combined.neurons <- FindVariableFeatures(object =  multiome.combined.neurons, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x =  multiome.combined.neurons)
multiome.combined.neurons <- ScaleData(object =  multiome.combined.neurons, features = all.genes)
multiome.combined.neurons <- RunPCA(object =  multiome.combined.neurons, features = VariableFeatures(object =  multiome.combined.neurons))
ElbowPlot(object =  multiome.combined.neurons, ndims = 50)
multiome.combined.neurons <- FindNeighbors(object =  multiome.combined.neurons, dims = 1:20)
multiome.combined.neurons <- FindClusters(object =  multiome.combined.neurons, resolution = 1.0) #2.0 before
multiome.combined.neurons <- RunUMAP(object =  multiome.combined.neurons, dims = 1:20, return.model=TRUE)
DimPlot(object =  multiome.combined.neurons, reduction = "umap", label=TRUE, pt.size = 0.25)
DefaultAssay(multiome.combined.neurons) <- "RNA"


#save(multiome.combined.neurons, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Rcodes/CombinedMultiomeData_250110.RData")


## Align to 2020 UMAP

aggregate.combined <- UpdateSeuratObject(aggregate.combined)
aggregate.combined <- RunUMAP(object = aggregate.combined, reduction = "pca", dims = 1:20, return.model=TRUE)
aggregate.combined$renamed_clusters <-Idents(aggregate.combined)
multiome_anchors <- FindTransferAnchors(reference = aggregate.combined, query = multiome.combined.neurons, dims = 1:30, reference.reduction = "pca")
multiome.combined.neurons <- MapQuery(anchorset = multiome_anchors, query = multiome.combined.neurons, reference = aggregate.combined, refdata = list(predicted_clusters = "renamed_clusters"), reference.reduction = "pca", reduction.model = "umap")


DimPlot(multiome.combined.neurons, group.by = "predicted.predicted_clusters",reduction = "ref.umap", label = T)


get_condition <- function(cell_name) {
  if (startsWith(cell_name, "NT_")) return("NT")
  if (startsWith(cell_name, "PBS_C_")) return("PBS_C")
  if (startsWith(cell_name, "OVA_C_")) return("OVA_C")
  if (startsWith(cell_name, "PBS_")) return("PBS")
  if (startsWith(cell_name, "OVA_")) return("OVA")
  return(NA)
}

multiome.combined.neurons$condition <- sapply(colnames(multiome.combined.neurons), get_condition)





DimPlot(multiome.combined.neurons, group.by = "condition",reduction = "ref.umap", cols = c("blue", "orange", "red", "purple", "green"), pt.size = 0.5)


save(multiome.combined.neurons, file="C:/Users/ekroj/Dropbox (MIT)/PrescottLab/Rcodes/CombinedMultiomeData_250110.RData")





