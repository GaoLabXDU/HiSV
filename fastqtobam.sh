fastq1="/home/li/NAS/Hi-C_Data/Cell_lines/NCI-H460/ENCFF136ZMZ.fastq"
fastq2="/home/li/NAS/Hi-C_Data/Cell_lines/NCI-H460/ENCFF401FJO.fastq"
ref="/media/li/Data/hg38/hg38.fa"
refrag="/media/li/Data/hg38/restriction_site/hg38_hindiii.bed"
outDir="/home/li/NAS/Hi-C_Data/Cell_lines/NCI-H460"
binsize=50000
coreN=16

# echo "I. 1. alignment fastq_1"
# bwa mem -t $coreN $ref $fastq1 > $outDir/sample_1.sam
echo "I. 2. alignment fastq_2"
bwa mem -t $coreN $ref $fastq2 > $outDir/sample_2.sam
echo "II. buliding bam file"
hicBuildMatrix --samFiles $outDir/sample_1.sam $outDir/sample_2.sam --binSize $binsize --restrictionSequence AAGCTT --danglingSequence AGCT --restrictionCutFile $refrag --threads 4 --inputBufferSize 100000 --outBam $outDir/sample.bam -o $outDir/sample_matrix.h5  --QCfolder $outDir/hicQC
