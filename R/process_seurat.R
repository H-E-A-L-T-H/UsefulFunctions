
#' Process Seurat objects
#'
#' This command creates a Seurat object, filters for doublets and dead cells
#' runs UMap, Tsne, and creates a clusterplot
#' @param dirlist takes a list of 10x sequencing directories to process
#' @param outlist takes a list of what you would like the output objects to be called
#' @keywords Seurat processing
#' @export
#' @return a list of Seurat objects that in which doublets and dead cells have been filtered
#' variable features have been identified, normalized by SCtransform,
#' @details you may need to modify the criteria for subsetting. If your elbow plot is far from
#' 18 you will need to change it most likely.
#' process_seurat()

process_seurat <- function(dirlist, outlist){
  lapply(dirlist, FUN = function(λ) {
    temp <- Seurat::Read10X(data.dir = λ)
    temp <- Seurat::CreateSeuratObject(counts = temp, outlist[dirlist == λ])
    temp$model <- outlist[λ]
    temp[["percent.mt"]] <- Seurat::PercentageFeatureSet(temp, pattern = "mt")
    Seurat::VlnPlot(temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    temp <- subset(temp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    temp <- Seurat::FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
    temp <- Seurat::SCTransform(temp, verbose = F)
    temp <- Seurat::RunPCA(temp, verbose = FALSE)
    Seurat::ElbowPlot(temp)

    temp <- Seurat::RunUMAP(temp, dims = 1:18, verbose = FALSE)
    temp <- Seurat::RunTSNE(temp, dims = 1:18, verbose = FALSE)
    temp <- Seurat::FindNeighbors(temp, dims = 1:18, verbose = FALSE)
    temp <- Seurat::FindClusters(temp, verbose = FALSE)
    assign(outlist[dirlist == λ],temp, envir = .GlobalEnv)
  })
}
