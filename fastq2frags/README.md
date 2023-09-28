## workflow

<div>
    <img src="../static/snakemake_dag.png" alt="dag" width="400"/>
</div>

The fastq2frags pipeline requires the barcoding schema described in [Zhang et al. 2023](https://www.nature.com/articles/s41586-023-05795-1#MOESM1) for Illumina short-read sequencing:
- read1: genomic sequence
- read2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence

1. Filter reads on read2 ligation linker sequences with [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/).
    1. reads with >3 mismatches in the ligation linker1 are removed from analysis
    2. reads with >3 mismatches in the ligation linker2 are removed from analysis

2. Perform sequence alignment with [Chromap](https://github.com/haowenz/chromap).

3. Convert the output file from the [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) to a fragments.tsv.gz file.
    
## setting up the environment
The ‘fastq2frags' pipeline can be run in an AWS EC2 instance or other cloud computing resource.  We developed the pipeline in an EC2 instance running Red Hat Enterprise Linux 9 (RHEL 9) with 32 cores and 64GB of RAM (ami-08e637cea2f053dfa).  The pipeline requires at least 1 TB free of free disk space.  See below for instructions for setting up an environment for fastq2frags in RHEL 9.

> For an example of running the pipeline in a high performance cluster (HPC) see [here](https://github.com/di-0579/Spatial_epigenome-transcriptome_co-sequencing/tree/main/Data_preprocessing/Spatial-ATAC-seq).

1. Install `git` and set up ssh [access](https://www.theodinproject.com/lessons/foundations-setting-up-git).
    ```
    sudo dnf update 
    sudo dnf install git
    ```

2. Install other dependencies. Java v-7 or greater is required to run bbmap. 
    ```
    sudo dnf install curl gcc gcc-c++ java-11-openjdk-devel make unzip zlib-devel
    ```

3. Install `bgzip`; in RHEL 9, `bgzip` in included as part of [htslib](http://www.htslib.org/download/).

4. Clone the [ATC_ATAC-seq](https://github.com/atlasxomics/ATX_ATAC-seq/) repository.
    ```
    git clone git@github.com:atlasxomics/ATX_ATAC-seq.git 
    ```

5. In the ATX_ATAC-seq/fastq2frags directory, download and install [chromap](https://github.com/haowenz/chromap#install).
    ```
    curl -L https://github.com/haowenz/chromap/archive/refs/heads/li_dev4.zip -o li_dev4.zip
    unzip li_dev4.zip
    mv chromap-li_dev4 chromap
    cd chromap && make
    ```

6. Download a reference genome to the fastq2frags directory and create an index for Chromap [index](https://github.com/haowenz/chromap#general).
For access to premade indices and reference genomes for human (GRCh38) and mouse (GRCm38), please contact your AtlasXomics support scientist or 
email (support@atlasxomics.com).

7. Download [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/) into the fastq2frags repository.
    ```
    curl -L https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download -o BBMap_39.01.tar.gz 
    tar -xvzf BBMap_39.01.tar.gz 
    ```

8. Install [mambaforge](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1-installing-mambaforge).
    ```
    curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
    bash Mambaforge-Linux-x86_64.sh
    ```

9. Install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) into a new conda enviroment.
    ```
    conda activate base
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    ```

When the environment has been correctly setup, the fastq2frags directory should have the following structure:

<div>
    <img src="../static/fastq2frags.png" alt="fastq2frags" width="400"/>
</div>

## running the pipeline
1. Upload fastq files into the fastqs directory in a subdirectory named with the sample name: ie. fastq2frags/fastq/SAMPLENAME/[SAMPLENAME_R1.fastq.gz, SAMPLENAME_R2.fastq.gz]

- AtlasXomics delivers fastq files via [latch.bio](https://latch.bio/); to facilitate transferring data from latch, you may want to install **latch cli**.
    ```
    python3 -m pip install latch
    ```
  Instructions for coping data from latch can be found [here](https://wiki.latch.bio/wiki/data/data-command-line).

2. Ensure the global variables in the Snakefile have the correct values.
* Check that `FASTQ_DIR` points to the top directory containing fastqs, with subdirectories for each sample containing R1 and R2 fastq.gz files.
* Check that `REF_DIR` points to the correct reference genome directory, containing a reference fasta (.fa or .fasta) and index (.index) file.
* Check that `BARCODE_FILE` points to the correct .txt.gz barcode 'whitelist' file.  Use bc50.txt.gz or bc96.txt.gz for 50-channel and 96-channel chips, respectively.
* Check that `CORES` contains the correct number of available cores for the pipeline.  Here, we set the value to 32.

4. Activate the snakemake environment and run the workflow.  The parameter -c can be set the a specificed number of cores (ie. 32); see snakemake [docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html#command-line-interface).
    ```
    conda activate snakemake
    snakemake -c
    ```

5. Once the workflow has completed, the fastq2frags directory should have the following structure.  The fragments.tsv.gz file can be used directly in downstream processing (see the [analysis](../analysis/) directory for a tutorial on downstream processing).

<div>
    <img src="../static/dirs.png" alt="dag" width="400"/>
</div>

## support
Questions? Comments?  Contact support@atlasxomics.com or post in the AtlasXomics [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).
