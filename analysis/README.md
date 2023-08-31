# Spatial analysis of spatial, single cell epigenomic data

This repository contains R Markdown tutorials for the analysis of spatial
ATAC-seq data generated via DBIT-seq.

These tutorials assume you have created fragments files from a epigenomic 
alignment and preprocessing pipeline (ie.[Chromap](https://www.nature.com/articles/s41467-021-26865-w),
[Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac)),
and 'spatial' folders via our [AtlasXBrowser](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/Overview.html)
software.  For our custom alignment and preprocessing workflow, see `fastq2frags/` in
this repository. 

This tutorial assumes that example data (fragments, spatial folders) is saved in
this directory with the following structure:

<div>
    <img src="./figures/tree.png" alt="dag" width="300"/>
</div>

To access example data, please contact your AtlasXomics support scientist or contact support@atlasxomics.com.

## dependiencies
* [ArchR v1.0.2](https://www.archrproject.com/)
    > ArchR is designed to run on UNIX-based platforms (Linux, Mac); it is not support for Windows.
* [Seurat v4.3.0](https://satijalab.org/seurat/)
* [macs2 v2.2.6](https://pypi.org/project/MACS2/2.2.6/)
