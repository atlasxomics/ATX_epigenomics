# Analysis of an epigenomic DBiT-seq experiment

This directory contains an R Markdown [tutorial](epigenomics_tutorial.Rmd) describing common tasks performed
in the analysis of spatial epigenomic data generated via DBiT-seq.  This analysis relies heavily on the R packages 
[ArchR](https://www.archrproject.com/) and [Seurat](https://satijalab.org/seurat/articles/install_v5#install-seurat-v4);
we recommend reading their vingettes as a companion to our tutorial.

This tutorial assumes you have created fragments files from a epigenomic alignment and preprocessing pipeline (ie.
[Chromap](https://www.nature.com/articles/s41467-021-26865-w)), and 'spatial' folders via our [AtlasXBrowser](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/Overview.html) software.  For our custom alignment and
preprocessing workflow, see `fastq2frags/` in this repository. 

This tutorial assumes that example data (fragments, spatial folders) is saved in this directory with the following
structure:

<div>
    <img src="./figures/tree.png" alt="dag" width="300"/>
</div>

To access example data, please contact your AtlasXomics support scientist or contact support@atlasxomics.com.

## Summary
In this tutorial, we recreeate the analysis of a spatial CUT&Tag experiment performed internally by 
AtlasXomics.  We processed three adjecnt fresh frozen mouse hippocampal brain sections (D1208, D1209, D1210)
according to  [Deng, 2022](https://www.science.org/doi/10.1126/science.abg7216), with an antibody for H3K27ac.
The three sections were treated as replicates.  We sequenced the NGS libraries on an Illumina NextSeq with a
target depth of 200 million reads, and transformed the resultant fastq files into epigenomic fragment files
via the [fastq2frag](../fastq2frags/) pipeline.

We start this tutorial by creating Arrow files and an ArchRProject from the fragment.tsv.gz files.
'Off-tissue' tixels (cells) are filtered from the ArchRProject project.  We then perform dimensionality
reduction and clustering with functions from ArchR.  Metadata and gene scores are extracted from ArchRProject
and combined with image data from the spatial folders to create SeuratObjects for each tissue.  The 
SeuartObjects are used to project data onto the tissue, allowing for spatial analysis of the data. Following
standard single-cell and epigenomic analyses workflows, we then perform peaks calling, and identify 
differentially expressed peaks, genes, and motifs.  Finally, we approxamate cell-typing by mapping the 
putative expression of marker genes onto the tissues. 

## Environment Setup
We recommend using RStudio for this tutorial; we've successsfully run the code with R v4.3.1+.  We recommend running
on a Unix-based machine (Mac, Linux) with at least 16GB of RAM and 4 CPUs.  Before starting the analysis, ensure you 
have the example data downloaded into this directory, as specified in the figure above.

This tutorial was written with the speicific versions of the following packages:
* [ArchR v1.0.2](https://www.archrproject.com/)
    > ArchR is designed to run on UNIX-based platforms (Linux, Mac); it is not supported for Windows.
* [Seurat v4.3.0](https://satijalab.org/seurat/articles/install_v5#install-seurat-v4)
* [SeuratObject v4.1.1](https://github.com/satijalab/seurat-object/releases/tag/v4.1.3)
* [macs2 v2.2.6](https://pypi.org/project/MACS2/2.2.6/)
    * We have found peak calling requires this specifc version of macs2, which relies on numpy v1.26.2.  You
    can install them with,
      ```
      pip install numpy==1.26.2
      pip install MACS2==2.2.6
      ```

Regreatablly, we have found these packages are susecptable to compatibility issues.  To facilitate code 
reproduciblability, we created an [renv](https://rstudio.github.io/renv/articles/renv.html) for this tutorial.

> Since this tutorial was written, SeuratV5 has been release.  This version is incompatible with the code in this
> tutorial.  We are working on an updated tutorial for SeuratV5.

### Instructions for renv

1. Clone this respoitory, open RStudio, and set your working directory to `ATX_epigenomics
/analysis/`. 
2. Install `renv`.
   ```
   install.packages('renv')
   ```
3. Activate `renv`.
   ```
   renv::activate()
   ```
5. Initiate setting up the project R environment; selet 'Y' to proceed with install and downloading the R packages
   into tutorial project.
   ```
    renv::restore()
   ```
   * Some R package may have system library dependencies (i.e. libgsl-dev) that need to be installed. If so, `renv::restore()`
     will break with an error message describing the missing dependencies.  Install the dependences and resume the enviroment
     setup by re-running `renv::restore()`.
6.    


