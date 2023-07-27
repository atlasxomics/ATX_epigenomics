## pipeline workflow

<div align="center">
    <img src="../static/snakemake_dag.png" alt="dag" width="400"/>
</div>

> This pipeline was built for short-read NGS data generated via Illumina platforms; the pipeline may require modification for NGS data from other platforms.

For this pipeline expects the following barcoding schema,
- read 1: genomic sequence
- read 2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence

### steps

1. filter reads on reads on read2 ligation linker sequences with [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
    1. reads with >3 mismatchs in the ligation linker1 removed from analysis
    2. reads with >3 mismatchs in the ligation linker2 removed from analysis

2. reformat read2 fastq to Cell Ranger ATAC format with BC_process.py
    - read2 -> read3, new-read2
        - read3: genome sequence
        - new-read2: barcodeA + barcodeB

3. sequence alignment and fragment file generation with Cell Ranger ATAC
    
    The following reference genomes were used,
    
    - Mouse reference (mm10)
        ```
        curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.2.0.tar.gz
        ```

    - Human reference (GRCh38):
        ```
        curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz
        ```
## setting up the enviroment
Due to the system requirements of Cell Ranger ATAC, the pipeline must be run in a Linux environment with at least 64GB RAM, an 8-core Intel or AMD processor, and 1TB of free disk space; we recommend running this pipeline in an AWS EC2 instance.

An example of running the pipeline in an HPC can be found [here](https://github.com/di-0579/Spatial_epigenome-transcriptome_co-sequencing/tree/main/Data_preprocessing/Spatial-ATAC-seq).

## running the pipeline
<div>
    <img src="../static/dirs.png" alt="dag" width="400"/>
</div>

1. Ensure cellranger-atac-cs/1.2.0/lib/python/barcodes/737K-cratac-v1.txt matches the barcode file in top directory.
2. Configure Snakefile.
3. Activate conda snakemake env and run `snakemake -j 32`

## dependiencies

* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
* [Biopython](https://biopython.org/docs/1.75/api/index.html)
* [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation). v1.2
* [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
