<div align="center">

![data](images/data.png)

# Spatial profiling of chromatin accessibility via DBIT-seq

</div> 

Deterministic Barcoding in Tissue for spatial omics sequencing (DBiT-seq) uncovers spatial biology by combining microfluidics and next- generation sequencing (NGS). Applications of the platform published in [Nature](https://www.nature.com/articles/s41586-022-05094-1), [Nature Biotech](https://www.nature.com/articles/s41587-023-01676-0), [Science](https://www.science.org/doi/10.1126/science.abg7216), and [Cell](https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420313908%3Fshowall%3Dtrue), demonstrate the platform’s comprehensive, unbiased profiling of the transcriptome, proteome and, for the first time, epigenome atthe cellular level.

This repository contains the  scripts necessary to analyze data from ATAC DBIT-seq experiments performed with the barcoding scheme from [Zhang et al. 2023](https://www.nature.com/articles/s41586-023-05795-1#MOESM1) (barcodes on read2).

<div align="center">
  
[AtlasXomics](https://www.atlasxomics.com) • [Docs](https://docs.atlasxomics.com) • [Discord](https://discord.com/channels/1004748539827597413/1004748540624511008)

</div> 

# Workflow
Workflows developed with the SDK feature:

  * Instant no-code interfaces for accessibility and publication
  * First class static typing
  * Containerization + versioning of every registered change
  * Reliable + scalable managed cloud infrastructure
  * Single line definition of arbitrary resource requirements (eg. CPU, GPU) for serverless execution

The Latch SDK is a framework to build workflows. A collection of existing and
maintained workflows for common biological assays can be found at [Latch
Verified](https://github.com/latch-verified).

### preprocessing

1. Raw Fastq data

    * Read 1: contains the spatial Barcode A and Barcode B

    * Read 2: contains the genome sequences

2. Reformat raw Fastq file to Cell Ranger ATAC format (10x Genomics)

    Reformatting raw data was implemented by BC_process.py in the Data_preprocessing folder.
   
    * Raw read 1 -> New Read 1 + New Read 2**

    * New Read 1: contains the genome sequences

    * New Read 2: contains the spatial Barcode A and Barcode B

    * **Raw read 2 -> New Read 3**

5. Sequence alignment and generation of fragments file

    The reformated data was processed using Cell Ranger ATAC v1.2 with following references:
    
    Mouse reference (mm10):
    ```
    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.2.0.tar.gz
    ```
    
    Human reference (GRCh38):
    ```
    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz
    ```
    
    **A preprocessing pipeline we developed using Snakemake workflow management system is in the Data_preprocessing folder. To run the pipeline, use the command:**
    ```
    sbatch Snakemake.sh
    ```

### analysis
  The data visualization were implemented with ArchR v1.0.1 and Seurat v3.2 package (Data_visualization folder).
  
  Brief descriptions of analysis scripts:
  
  **metadata_files_for_Seurat_spatial.ipynb**: Generate metadata files that were compatible with Seurat workflow for spatial datasets.
  
  **archR.R**: Data normalization and dimensionality reduction, identifying the marker genes, peak calling, deviatons enrichment anaylsis, bulk sample projection, and pseudotime analysis.
  
  **spatial_data_visualization.R**: Visualize spatially resolved data on tissue sections.
  
  **GO_enrichment_analysis.R**: GO enrichment analysis for marker genes.
  
  **integrative_data_analysis.R**: Integrative data analysis with scRNA-seq reference datasets.

# Citations

Deng, Y., Bartosovic, M., Ma, S. et al. Spatial profiling of chromatin accessibility in mouse and human tissues. Nature 609, 375–383 (2022). https://doi.org/10.1038/s41586-022-05094-1

Liu, Y., DiStasio, M., Su, G. et al. High-plex protein and whole transcriptome co-mapping at cellular resolution with spatial CITE-seq. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01676-0

Yanxiang Deng et al. ,Spatial-CUT&Tag: Spatially resolved chromatin modification profiling at the cellular level.Science375,681-686(2022). DOI:10.1126/science.abg7216

Liu et al. High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue Cell. 183, 1665–1681 (2020). https://doi.org/10.1016/j.cell.2020.10.026




