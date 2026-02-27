## workflow

The fastq2frags pipeline requires the barcoding schema described in [Zhang et al. 2023](https://www.nature.com/articles/s41586-023-05795-1#MOESM1) for Illumina short-read sequencing:
- read1: genomic sequence
- read2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence

1. Filter reads on read2 ligation linker sequences with [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/).
    1. reads with >3 mismatches in the ligation linker1 are removed from analysis
    2. reads with >3 mismatches in the ligation linker2 are removed from analysis

2. Perform sequence alignment with [Chromap](https://github.com/haowenz/chromap).

3. Convert the output file from the [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) to a fragments.tsv.gz file.

<div align="center">
    <img src="../static/snakemake_dag.png" alt="dag" width="175"/>
</div>
    
## setting up the environment
The `fastq2frags` pipeline can be run locally or in an AWS EC2 instance or other cloud computing resource. We developed the pipeline in an EC2 instance running Red Hat Enterprise Linux 9 (RHEL 9) with 32 cores and 64GB of RAM (ami-08e637cea2f053dfa). We recommend at least 1 TB free of disk space for full experiments, but less can be used for testing. See below for instructions on setting up an environment for fastq2frags in RHEL 9.

**relevant package versions**
- Chromap: li_dev4
- BBMap: BBMap_39.01

> For an example of running the pipeline in a high-performance cluster (HPC) see [here](https://github.com/di-0579/Spatial_epigenome-transcriptome_co-sequencing/tree/main/Data_preprocessing/Spatial-ATAC-seq). AtlasXomics does not support this repository.

1. Ensure `git` is installed.
    ```
    sudo dnf update 
    sudo dnf install git
    ```
    or, for Ubuntu,
    ```
    sudo apt-get update -y
    sudo apt-get install -y git
    ```

2. Install other dependencies. Java v-7 or greater is required to run bbmap.
    ```
    sudo dnf install curl gcc gcc-c++ java-11-openjdk-devel make unzip zlib-devel
    ```
    or
    ```
    sudo apt-get install -y curl default-jdk build-essential make unzip zlib1g-dev wget
    ```

3. Install `bgzip`. In RHEL 9, `bgzip` is included as part of htslib; see download instructions [here](http://www.htslib.org/download/). For Ubuntu, it is part of the package `tabix`.
    ```
    sudo apt-get install tabix
    ```

4. Clone the [ATX_epigenomics](https://github.com/atlasxomics/ATX_epigenomics) repository.
    ```
    git clone https://github.com/atlasxomics/ATX_epigenomics.git
    ```
    Once the repository has been cloned, navigate to `ATX_epigenomics/fastq2frags`. 
    ```
    cd ATX_epigenomics/fastq2frags/
    ```

5. Download and install [chromap](https://github.com/haowenz/chromap#install) into the `fastq2frags` directory.
    ```
    curl -L https://github.com/haowenz/chromap/archive/refs/heads/li_dev4.zip -o li_dev4.zip
    unzip li_dev4.zip
    mv chromap-li_dev4 chromap
    cd chromap && make
    ```
    Once `chromap` has been installed, navigate back to `ATX_epigenomics/fastq2frags`. 
    ```
    cd ..
    ```

6. Download a reference genome to the fastq2frags directory. Premade indices and reference genomes for human (GRCh38), mouse (GRCm38), and rat (Rnor6) can be downloaded via the links below:
- GRCh38:
    - GUI: https://console.latch.bio/s/17114493928624682
    - URL: https://latch-public.s3.amazonaws.com/test-data/13502/chromap_references/GRCh38.tar.gz
- GRCm38
    - GUI: https://console.latch.bio/s/17114488994250799
    - URL: https://latch-public.s3.amazonaws.com/test-data/13502/chromap_references/GRCm38.tar.gz
- Rnor6
    - GUI: https://console.latch.bio/s/17114503950478690
    - URL: https://latch-public.s3.amazonaws.com/test-data/13502/chromap_references/Rnor6.tar.gz

    Compressed tar references can be downloaded with `curl` or `wget` and extracted:
    ```
    wget https://latch-public.s3.amazonaws.com/test-data/13502/chromap_references/GRCm38.tar.gz
    tar -xzvf GRCm38.tar.gz
    ```

    Set `REF_DIR` in `Snakefile` to the extracted reference directory name (for example `REF_DIR = 'GRCm38'` if you extracted `GRCm38.tar.gz`).

    See Chromap documentation for how to prepare a custom genome [index](https://github.com/haowenz/chromap#general).

7. Download [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/) into the fastq2frags repository.
    ```
    curl -L https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download -o BBMap_39.01.tar.gz 
    tar -xvzf BBMap_39.01.tar.gz 
    ```

8. Install [mambaforge](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1-installing-mambaforge).
    ```
    curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh
    ```

9. Install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) into a new conda environment.
    ```
    conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
    ```
    The snakemake environment can be activated with:
    ```
    conda activate snakemake
    ```

10. Make a directory called `fastqs`. Within `fastqs`, make a subfolder named with your sample name. If analyzing more than one sample, make a directory for each. Move or download FASTQ files for the sample into the corresponding sample subfolder.
    ```
    mkdir fastqs
    mkdir fastqs/test_dataset_50m
    ```

When the environment has been correctly setup, the fastq2frags directory should have the following structure:

```text
fastq2frags/
├── bbmap/
├── chromap/
├── bc50.txt.gz
├── bc96.txt.gz
├── Snakefile
├── <REF_DIR>/
│   ├── *.fa
│   └── *.index
└── fastqs/
    └── SAMPLENAME/
        ├── SAMPLENAME_R1.fastq.gz
        └── SAMPLENAME_R2.fastq.gz
```

## running the pipeline
1. Upload FASTQ files into the `fastqs` directory in a subdirectory named with the sample name:
   `fastq2frags/fastqs/SAMPLENAME/[SAMPLENAME_R1.fastq.gz, SAMPLENAME_R2.fastq.gz]`

- AtlasXomics delivers fastq files via [latch.bio](https://latch.bio/); to facilitate transferring data from latch, you may want to install **latch cli**.
    ```
    python3 -m pip install latch
    ```
  Instructions for copying data from latch can be found [here](https://wiki.latch.bio/wiki/data/data-command-line).

2. Ensure the global variables in the Snakefile have the correct values.
* Check that `FASTQ_DIR` points to the top directory containing fastqs, with subdirectories for each sample containing R1 and R2 fastq.gz files.
* Check that `REF_DIR` points to the correct reference genome directory, containing a reference fasta (.fa or .fasta) and index (.index) file.
* Check that `BARCODE_FILE` points to the correct .txt.gz barcode 'whitelist' file.  Use bc50.txt.gz or bc96.txt.gz for 50-channel and 96-channel chips, respectively.
* Check that `CORES` contains the correct number of available cores for the pipeline.  Here, we set the value to 32.
* Check that `BBDUK_XMX` sets a suitable Java heap size for bbduk (for example, `BBDUK_XMX = '24g'`).

3. Activate the snakemake environment and run the workflow. The parameter `-c` can be set to a specified number of cores (e.g., 8); see snakemake [docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html#command-line-interface).
    ```
    conda activate snakemake
    snakemake -c 8
    ```

4. Once the workflow has completed, outputs are organized as `outputs/{sample}/`. The `*_fragments.tsv.gz` file is the final output for downstream processing (see the [analysis](../analysis/) directory for a tutorial). Intermediate files listed below are workflow artifacts.

    ```text
    fastq2frags/
    ├── fastqs/
    │   └── SAMPLENAME/
    │       ├── SAMPLENAME_R1.fastq.gz
    │       └── SAMPLENAME_R2.fastq.gz
    └── outputs/
        └── SAMPLENAME/
            ├── SAMPLENAME_linker1_R1.fastq.gz
            ├── SAMPLENAME_linker1_R2.fastq.gz
            ├── SAMPLENAME_linker2_R1.fastq.gz
            ├── SAMPLENAME_linker2_R2.fastq.gz
            ├── SAMPLENAME_stats.linker1.txt
            ├── SAMPLENAME_stats.linker2.txt
            ├── SAMPLENAME_aln.bed
            ├── SAMPLENAME_temp.bed
            ├── SAMPLENAME_fragments.tsv
            └── SAMPLENAME_fragments.tsv.gz
    ```

## support
Questions? Comments?  Contact support@atlasxomics.com.
