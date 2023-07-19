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

**preprocessing** contains a scATAC-seq analysis pipeline built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home); raw fastq files are uploaded directly to this directory for processing.  The pipeline performes the following steps,

- reads are filtered on the correctness of ligation linker sequences
- read2 is split into two read files ('new' read2, read3)
- Cell Ranger ATAC is run on the processed reads.

Cell Ranger ATAC outputs fragment files that are used for downstream analysis.  Due to the [system requirements](https://support.10xgenomics.com/single-cell-atac/software/overview/system-requirements) of Cell Ranger ATAC, the pipeline must be run in a Linux enviroment with at least 64GB RAM, an 8-core Intel or AMD processor, and 1TB of free disk space.

**preprocessing** contains R scripts and vignettes for further analysis of ATAC DBIT-seq experiments.  This scripts rely on the [ArchR](https://www.nature.com/articles/s41588-021-00790-6) and Seurat packages.

# Citations

Deng, Y., Bartosovic, M., Ma, S. et al. Spatial profiling of chromatin accessibility in mouse and human tissues. Nature 609, 375–383 (2022). https://doi.org/10.1038/s41586-022-05094-1

Liu, Y., DiStasio, M., Su, G. et al. High-plex protein and whole transcriptome co-mapping at cellular resolution with spatial CITE-seq. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01676-0

Yanxiang Deng et al. ,Spatial-CUT&Tag: Spatially resolved chromatin modification profiling at the cellular level.Science375,681-686(2022). DOI:10.1126/science.abg7216

Liu et al. High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue Cell. 183, 1665–1681 (2020). https://doi.org/10.1016/j.cell.2020.10.026

Granja, J.M., Corces, M.R., Pierce, S.E. et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nat Genet 53, 403–411 (2021). https://doi.org/10.1038/s41588-021-00790-6






