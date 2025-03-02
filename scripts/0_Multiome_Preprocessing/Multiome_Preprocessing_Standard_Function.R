library(Seurat)
library(Signac)
library(GenomicFeatures)
library(EnsDb.Mmusculus.v79)
library(Rsamtools)

# Define a function that encapsulates creating the combined RNA+ATAC Seurat object
createRNA_ATAC_Seurat <- function(feature_matrix_dir, 
                                  fragments_path,
                                  genome_name = "mm10",
                                  annotation_db = EnsDb.Mmusculus.v79,
                                  pattern_mt = "^mt-",
                                  min_cells_atac = 10) {
  # 1. Read in data from 10x directory (contains both RNA and ATAC).
  data <- Read10X(data.dir = feature_matrix_dir)
  
  # 2. Extract RNA and ATAC assays
  rna_counts  <- data$`Gene Expression`
  atac_counts <- data$Peaks
  
  # 3. Create Seurat object with RNA counts
  seurat_obj <- CreateSeuratObject(counts = rna_counts)
  
  # 4. Calculate percent.mt (optional but typically useful)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern_mt)
  
  # 5. Prepare ATAC data: filter out non-standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts   <- atac_counts[as.vector(grange.use), ]
  
  # 6. Prepare genomic annotations
  annotations <- GetGRangesFromEnsDb(ensdb = annotation_db)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- genome_name
  
  # 7. Create ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts     = atac_counts,
    sep        = c(":", "-"),
    genome     = genome_name,
    fragments  = fragments_path,
    min.cells  = min_cells_atac,
    annotation = annotations
  )
  
  # 8. Add the ATAC assay to the Seurat object
  seurat_obj[["ATAC"]] <- chrom_assay
  
  # 9. Return the combined object
  return(seurat_obj)
}

