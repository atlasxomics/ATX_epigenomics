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

spatial_plot <- function(seurat_object, name) {
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

feature_plot <- function(seurat_obj, feature, name) {
  SpatialFeaturePlot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1
  ) +
    ggtitle(paste0(feature, " : ", name)) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 15)
    )
}
