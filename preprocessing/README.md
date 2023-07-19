# Spatial ATAC-seq
A snakemake pipeline to process spatial ATAC-seq raw data

## Work flow of the pipeline

![](./snakemake_dag.png)

## Dependiencies

* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
* [Biopython](https://biopython.org/docs/1.75/api/index.html)
* [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation). v1.2
* [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)

Next Generation Sequencing (NGS) was performed using the Illumina NextSeq 2000 or NovaSeq 6000 (150:8:150). 

### 1. Raw Fastq data

Read 1: contains the genomic sequence

Read 2: contains the spatial Barcode A and Barcode B and genomic sequence

### 2. Reformat raw Fastq file to Cell Ranger ATAC format (10x Genomics)

**Raw read 2 -> New Read 3 + New Read 2**

- New Read 3: contains the genome sequences

- New Read 2: contains the spatial Barcode A and Barcode B

**Raw read 2 -> New Read 3**

Reformatting raw data was implemented by BC_process.py in the Data_preprocessing folder.

### 3. Sequence alignment and generation of fragments file

The reformated data was processed using Cell Ranger ATAC v1.2 with following references:

Mouse reference (mm10):
```
curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.2.0.tar.gz
```

Human reference (GRCh38):
```
curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz
```

**A preprocessing pipeline we developed using Snakemake workflow management system is in the Data_preprocessing folder.**

## Run the pipeline
1. Ensure cellranger-atac-cs/1.2.0/lib/python/barcodes/737K-cratac-v1.txt matches the barcode file in top directory.
2. Configure Snakefile
3. To run the pipeline, activate conda snakemake env and run the command:
```
snakemake -j 32
```

