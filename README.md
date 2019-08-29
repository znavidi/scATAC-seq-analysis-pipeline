# ATACseq_analysis

Here are the pipelines are being used for analyzing single cell ATAC-seq samples with snaptools and cellranger-atac tools.

# snaptools Pipeline:

Generating snap file for a data sample could be started from raw fastq files or a XXX_fragments.tsv file(which might exist for some samples like 10X data) as input, which will be included here.



## What is snap file:

snap (Single-Nucleus Accessibility Profiles) file is a hierarchically structured hdf5 file that is specially designed for storing single nucleus ATAC-seq datasets. A snap file (version 4) contains the following sessions: header (HD), cell-by-bin accessibility matrix (AM), cell-by-peak matrix (PM), cell-by-gene matrix (GM), barcode (BD) and fragment (FM).

HD session contains snap-file version, created date, alignment and reference genome information.
BD session contains all unique barcodes and corresponding meta data.
AM session contains cell-by-bin matrices of different resolutions (or bin sizes).
PM session contains cell-by-peak count matrix. PM session contains cell-by-gene count matrix.
FM session contains all usable fragments for each cell. Fragments are indexed for fast search. Detailed information about snap file can be found here.

# Analyzing 10X data with snaptools:

### Step 1. Barcode demultiplexing (This step commands is for 10X Genomics data. Other datasets could use other scripts for demultiplexing)
--snaptools provide a module dex-fastq to integrate the 10X barcode into the read name (run this tool on R1 and R3 for all library files)--

snaptools dex-fastq --input-fastq=atac_v1_E18_brain_cryo_5k_S1_L003_R3_001.fastq.gz --output-fastq=atac_v1_E18_brain_cryo_5k_S1_L003_R3_001.dex.fastq.gz --index-fastq-list atac_v1_E18_brain_cryo_5k_S1_L003_R2_001.fastq.gz

combine these two library

cat Library1_1_L001_R1_001.fastq.gz Library1_2_L001_R1_001.fastq.gz > Library1_L001_R1_001.fastq.gz

cat Library1_1_L001_R3_001.fastq.gz Library1_2_L001_R3_001.fastq.gz > Library1_L001_R3_001.fastq.gz

### Step 2. Index reference genome (snaptools)
snaptools index-genome  \
	--input-fasta=hg38.fa  \
	--output-prefix=hg38  \
	--aligner=bwa  \
	--path-to-aligner=/opt/biotools/bwa/bin/  \
	--num-threads=5

--run the rest of the pipeline using Library1_S1_L001_R1_001.fastq.dex.gz and Library1_S1_L001_R3_001.fastq.dex.gz--


### Step 3. Alignment (snaptools)
snaptools align-paired-end --input-reference=hg38.fa --input-fastq1=atac_v1_hgmm_500_S1_R1_001.dex.fastq.gz --input-fastq2=atac_v1_hgmm_500_S1_R3_001.dex.fastq.gz --output-bam=atac_v1_hgmm_500_S1.bam --aligner=bwa --read-fastq-command=zcat --min-cov=0 --num-threads=5 --if-sort=True --tmp-folder=./ --overwrite=TRUE


### Step 4. Pre-processing (snaptools).
This step generates snap file from aligned bam file:

snaptools snap-pre --input-file=atac_v1_hgmm_500_S1.bam --output-snap=atac_v1_hgmm_500_S1.snap --genome-name=hg38 --genome-size=hg38.chrom.sizes --min-mapq=30 --min-flen=0 --max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=FALSE --overwrite=True --min-cov=100 --verbose=True

Note: --keep-single argument must be TRUE if the data is single end and FALSE if the data is paired end!


### Step 5. Cell-by-bin matrix (snaptools)
snaptools snap-add-bmat --snap-file=atac_v1_hgmm_500_S1.snap --bin-size-list 1000 2000 5000 10000 --verbose=True

### Step 6. Analyzing snap file with SnapATAC R packages:
--analyzing for datasets without label--

Rscript --vanilla ~/projects/def-wanglab/znavidi/code/snapATAC_analysis_real.R ~/projects/def-wanglab/ATAC-seq-data/scATAC-seq/atac_v1_pbmc_10k/ atac_v1_pbmc_10k_S1 5000 http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz 20 atac_v1_pbmc_10k_singlecell.csv

--analyzing with datasets with label--

Rscript --vanilla ~/projects/def-wanglab/znavidi/code/snapATAC_analysis_sim.R ~/projects/def-wanglab/ATAC-seq-data/scATAC-seq/GSE74310/ GSE74310 5000 http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz 20



## Cell Ranger:

There is a complete doumentation of Cell Ranger tool: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac

### Installation
For installation follow instructions from this page:

https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation

I downloaded Cell Ranger ATAC - 1.1.0 (April 16, 2019) version from https://support.10xgenomics.com/single-cell-atac/software/downloads/latest and mm10 cell ranger reference files (stored in "ref" folder in "data" folder)

tar -xzvf cellranger-atac-1.1.0.tar.gz

export PATH=/cluster/home/znavidig/tool/cell_ranger/cellranger-atac-1.1.0:$PATH

### cell ranger count

manual: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count

~/tool/cell_ranger/cellranger-atac-1.1.0/cellranger-atac count --id=atac_v1_adult_brain_fresh_5k --reference=/cluster/projects/bwanggroup/single_cell/ATAC_data/ref/mm10/refdata-cellranger-atac-mm10-1.1.0 --fastqs=/cluster/projects/bwanggroup/single_cell/ATAC_data/scATACseq/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fastqs


