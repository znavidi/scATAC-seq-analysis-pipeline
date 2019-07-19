# ATACseq_analysis


snap atac analysis:

# Step 1. Barcode demultiplexing (This step commands is for 10X Genomics data. Other datasets could use other scripts for demultiplexing)
# snaptools provide a module dex-fastq to integrate the 10X barcode into the read name (run this tool on R1 and R3 for all library files)
snaptools dex-fastq --input-fastq=atac_v1_E18_brain_cryo_5k_S1_L003_R3_001.fastq.gz --output-fastq=atac_v1_E18_brain_cryo_5k_S1_L003_R3_001.dex.fastq.gz --index-fastq-list atac_v1_E18_brain_cryo_5k_S1_L003_R2_001.fastq.gz

# combine these two library
cat Library1_1_L001_R1_001.fastq.gz Library1_2_L001_R1_001.fastq.gz > Library1_L001_R1_001.fastq.gz
cat Library1_1_L001_R3_001.fastq.gz Library1_2_L001_R3_001.fastq.gz > Library1_L001_R3_001.fastq.gz

# Step 2. Index reference genome (snaptools)
snaptools index-genome  \
	--input-fasta=mm10.fa  \
	--output-prefix=mm10  \
	--aligner=bwa  \
	--path-to-aligner=/opt/biotools/bwa/bin/  \
	--num-threads=5

# run the rest of the pipeline using Library1_S1_L001_R1_001.fastq.dex.gz and Library1_S1_L001_R3_001.fastq.dex.gz
# Step 3. Alignment (snaptools)
snaptools align-paired-end --input-reference=../../ref/hg19/hg19.fa --input-fastq1=atac_v1_hgmm_500_S1_R1_001.dex.fastq.gz --input-fastq2=atac_v1_hgmm_500_S1_R3_001.dex.fastq.gz --output-bam=atac_v1_hgmm_500_S1.bam --aligner=bwa --read-fastq-command=zcat --min-cov=0 --num-threads=5 --if-sort=True --tmp-folder=./ --overwrite=TRUE

# Step 4. Pre-processing (snaptools).
snaptools snap-pre --input-file=GSE74310_merged.nsorted.bam --output-snap=../GSE74310.snap --genome-name=hg38 --genome-size=../../../ref/hg38/hg38.chrom.sizes --min-mapq=30 --min-flen=0 --max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=FALSE --overwrite=True --min-cov=100 --verbose=True

# Step 5. Cell-by-bin matrix (snaptools)
snaptools snap-add-bmat --snap-file=atac_v1_pbmc_5k_S1_001.snap --bin-size-list 1000 5000 10000 --verbose=True

# Step 6. Analyzing snap file with SnapATAC R packages:
# analyzing for datasets without label
Rscript --vanilla ~/projects/def-wanglab/znavidi/code/snapATAC_analysis_real.R ~/projects/def-wanglab/ATAC-seq-data/scATAC-seq/atac_v1_pbmc_10k/ atac_v1_pbmc_10k_S1 5000 http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz 20 atac_v1_pbmc_10k_singlecell.csv

# analyzing with datasets with label
Rscript --vanilla ~/projects/def-wanglab/znavidi/code/snapATAC_analysis_sim.R ~/projects/def-wanglab/ATAC-seq-data/scATAC-seq/GSE74310/ GSE74310 5000 http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz 20
