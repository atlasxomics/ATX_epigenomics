library(ArchR)
library(ggplot2)
library(harmony)
library(patchwork)
library(Seurat)

build_atlas_seurat_object <- function(
    run_id,
    matrix,
    metadata,
    spatial_path) {
  # Prepare and combine gene matrix, metadata, and image for SeuratObject
  # for runs within a project.
  image <- Read10X_Image(
    image.dir = spatial_path,
    filter.matrix = TRUE
  )
  metadata <- metadata[metadata$Sample == run_id, ]

  matrix <- matrix[, c(grep(pattern = run_id, colnames(matrix)))]
  matrix@Dimnames[[2]] <- metadata@rownames
  matrix <- CreateAssayObject(matrix)

  object <- CreateSeuratObject(
    counts = matrix,
    assay  = "Spatial",
    meta.data = as.data.frame(metadata)
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- "Spatial"
  object[["slice1"]] <- image
  return(object)
}

feature_plot <- function(seurat_obj, feature, name) {
  # Wrapper of Seurat's SpatialFeaturePlot with specific aesthetics

  SpatialFeaturePlot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1) +
      ggtitle(name) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))
}

spatial_plot <- function(seurat_object, name) {
  # Wrapper of Seurat's SpatialDimPlot with specific aesthetics

  clusters <- sort(unique(seurat_object$Clusters))
  colors <- ArchRPalettes$stallion2[seq_len(length(clusters))]
  names(colors) <- clusters
  SpatialDimPlot(
    seurat_object,
    group.by = "Clusters",
    pt.size.factor = 1,
    cols = colors,
    stroke = 0) +
    ggtitle(name) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 10),
      legend.position = "bottom")
}

geneset_plot <- function(seurat_obj, marker_genes, name, title) {
  # Return a Seurat SpatialFeaturePlot of the average expression score for a
  # set of genes.

  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = marker_genes,
    name = name
  )
  plot <- SpatialFeaturePlot(
    object = seurat_obj,
    pt = 1,
    features = paste0(name, 1)) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
}
