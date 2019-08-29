# ATACseq_analysis

Here are the pipelines are being used for analyzing single cell ATAC-seq samples with snaptools and cellranger-atac tools.


#### snaptools Pipeline:

Generating snap file for a data sample could be started from raw fastq files or a XXX_fragments.tsv file(which might exist for some samples like 10X data) as input, which will be included here.


##### What is snap file:

snap (Single-Nucleus Accessibility Profiles) file is a hierarchically structured hdf5 file that is specially designed for storing single nucleus ATAC-seq datasets. A snap file (version 4) contains the following sessions: header (HD), cell-by-bin accessibility matrix (AM), cell-by-peak matrix (PM), cell-by-gene matrix (GM), barcode (BD) and fragment (FM).

HD session contains snap-file version, created date, alignment and reference genome information.
BD session contains all unique barcodes and corresponding meta data.
AM session contains cell-by-bin matrices of different resolutions (or bin sizes).
PM session contains cell-by-peak count matrix. PM session contains cell-by-gene count matrix.
FM session contains all usable fragments for each cell. Fragments are indexed for fast search. Detailed information about snap file can be found here.


#### Analyzing 10X data with snaptools:

The pipline are run on https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_E18_brain_fresh_5k dataset from 10X as an example.

##### Step 1. Barcode demultiplexing (This step commands is for 10X Genomics data. Other datasets could use other scripts for demultiplexing)
snaptools provide a module dex-fastq to integrate the 10X barcode into the read name (run this tool on R1 and R3 for all library files).

	snaptools dex-fastq --input-fastq=atac_v1_E18_brain_cryo_5k_S1_L001_R1_001.fastq.gz --output-fastq=atac_v1_E18_brain_cryo_5k_S1_L001_R1_001.dex.fastq.gz --index-fastq-list atac_v1_E18_brain_cryo_5k_S1_L001_R2_001.fastq.gz

	snaptools dex-fastq --input-fastq=atac_v1_E18_brain_cryo_5k_S1_L001_R3_001.fastq.gz --output-fastq=atac_v1_E18_brain_cryo_5k_S1_L001_R3_001.dex.fastq.gz --index-fastq-list atac_v1_E18_brain_cryo_5k_S1_L001_R2_001.fastq.gz

combine these two libraries.

	cat atac_v1_E18_brain_cryo_5k_S1_L001_R1_001.dex.fastq.gz atac_v1_E18_brain_cryo_5k_S1_L002_R1_001.dex.fastq.gz > atac_v1_E18_brain_cryo_5k_R1.dex.fastq.gz

	cat atac_v1_E18_brain_cryo_5k_S1_L001_R3_001.dex.fastq.gz atac_v1_E18_brain_cryo_5k_S1_L002_R3_001.dex.fastq.gz > atac_v1_E18_brain_cryo_5k_R3.dex.fastq.gz

run the rest of the pipeline using atac_v1_E18_brain_cryo_5k_R1.dex.fastq.gz and atac_v1_E18_brain_cryo_5k_R3.dex.fastq.gz.


##### Step 2. Index reference genome (snaptools)
	snaptools index-genome  \
		--input-fasta=hg38.fa  \
		--output-prefix=hg38  \
		--aligner=bwa  \
		--path-to-aligner=/opt/biotools/bwa/bin/  \
		--num-threads=5


##### Step 3. Alignment (snaptools)
	snaptools align-paired-end --input-reference=hg38.fa --input-fastq1=atac_v1_E18_brain_cryo_5k_R1.dex.fastq.gz --input-fastq2=atac_v1_E18_brain_cryo_5k_R3.dex.fastq.gz --output-bam=atac_v1_E18_brain_cryo_5k.bam --aligner=bwa --read-fastq-command=zcat --min-cov=0 --num-threads=5 --if-sort=True --tmp-folder=./ --overwrite=TRUE


##### Step 4. Pre-processing (snaptools).
This step generates snap file from aligned bam file:

	snaptools snap-pre --input-file=atac_v1_E18_brain_cryo_5k.bam --output-snap=atac_v1_E18_brain_cryo_5k.snap --genome-name=hg38 --genome-size=hg38.chrom.sizes --min-mapq=30 --min-flen=0 --max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=FALSE --overwrite=True --min-cov=100 --verbose=True

Note: --keep-single argument must be TRUE if the data is single end and FALSE if the data is paired end!


##### Step 5. Cell-by-bin matrix (snaptools)
	snaptools snap-add-bmat --snap-file=atac_v1_E18_brain_cryo_5k.snap --bin-size-list 1000 2000 5000 10000 --verbose=True


##### Step 6. Analyzing snap file with SnapATAC R packages:
Analyzing for datasets without label:

	Rscript --vanilla ~/projects/def-wanglab/znavidi/code/snapATAC_analysis_real.R \

	  atac_v1_E18_brain_cryo_5k.bam/ \ (working directory) 

	  atac_v1_E18_brain_cryo_5k \ (name of snap file without .snap extension)

	  5000 \ (the bin size to analyze)

	  http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz \ (the black list of mous and http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz for human) 

	  20 \ (number of principals component to consider)

	  atac_v1_E18_brain_cryo_5k_singlecell.csv (Per Barcode metrics (CSV), could be downloaded from 10X site)

analyzing with datasets with label:

	Rscript --vanilla ~/projects/def-wanglab/znavidi/code/snapATAC_analysis_sim.R \

	GSE74310/ \ (working directory)

	GSE74310 \ (name of snap file without .snap extension)]

	5000 \ (the bin size to analyze)

	http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz \ (the black list address) 

	20 (number of principal components to consider)


#### Analyzing GEO data with snaptools:

GEO datasets that we have analayzed don't have cellranger canonical structured fastq or bam files and in order to be able to analyze them with snaptools we should add Barcode demultiplexing step ourself. So, we align the sampels to reference genome and then add 16 length barcode to the beginning of their sam files. Then convert them sam files to bam files using samtools and merge all samples to one bam file. The snap file will be created from this name sorted bam file.


##### Step 1. Download GEO datasets:


##### Step 2. Align fastq files:


##### Step 3. Generate and add unique barcodes for each single cell sample to sam files:


##### Step 4: Convert sam files to bam files:


##### Step 5: Merge bam files to XXX-merged.bam


##### Step 6: Sort merged bam file by name:


##### Step 7: Pre-processing (snaptools): 


##### Step 8. Cell-by-bin matrix (snaptools)


## Cell Ranger:

There is a complete doumentation of Cell Ranger tool: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac

### Installation
For installation follow instructions from this page:

https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation

Cell Ranger ATAC - 1.1.0 (April 16, 2019) version is downloaded from https://support.10xgenomics.com/single-cell-atac/software/downloads/latest and hg38 cell ranger reference files are stored in projects/def-wanglab/ATAC-seq-data/ref/refdata-cellranger-atac-GRCh38-1.1.0/ address.

tar -xzvf cellranger-atac-1.1.0.tar.gz

export PATH=/cluster/home/znavidig/tool/cell_ranger/cellranger-atac-1.1.0:$PATH

### Run cell ranger count

manual: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count

~/tool/cell_ranger/cellranger-atac-1.1.0/cellranger-atac count --id=atac_v1_hgmm_500_S1 
  / --reference=/cluster/projects/bwanggroup/single_cell/ATAC_data/ref/mm10/refdata-cellranger-atac-mm10-1.1.0 
  / --fastqs=atac_v1_adult_brain_fresh_5k_fastqs


