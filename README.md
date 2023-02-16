# HiSV
# Overview
HiSV is a computational pipeline for structural variation detection from Hi-C data

Structural variations (SVs) play an essential role in the evolution of human genomes and are associated with cancer genetics and rare disease.

High-throughput chromosome capture (Hi-C) is a new technology that breaks through the limitations of the short-read-based whole-genome sequenc-ing (WGS) method and can be used to detect large-scale SVs. However, the identification of SVs from Hi-C data is still in the early stages with only a few methods available. Especially, no algorithm has been developed that can detect major types of SVs without control samples.

Here, we described HiSV (Hi-C for Structural Variation), a control-free computational pipeline-based salient object detection model to identify large-scale SVs from a Hi-C sample. 

# Required Packages
The following software must be installed on your machine:

python: tested with version 3.8

## Python dependencies:
* numpy (tested 1.19.4)
* pysam (tested 0.16.0.1)
* pandas (tested 1.1.1)
* prox_tv (tested 3.3.0)

You can install the above package using the following command:
```
pip install numpy pysam pandas prox_tv
```

# Installation
Download HiSV by
```
git clone https://github.com/gaolabXDU/HiSV

```

# Usage
HiSV requires two input files, a tumor Hi-C data and a reference length file. HISV supports five types of Hi-C data as input including matrix file, bam file, .cool file, .hic file and fastq file. When Hi-C data is a bam file, .cool file, or hic file, you can use the script ``hiccovert.py`` to generate matrix file. When Hi-C data are fastq files, you can use the shell script ``fastqtobam.sh`` to generate bam file. This step requires the installation of BWA-MEM(https://github.com/lh3/bwa) and HiCExplorer(https://github.com/deeptools/HiCExplorer). HiSV results consist of two files including intra-chromosomal SV events and inter-chromosomal translocation events. 

## Ruuning command
```
Usage: HiSV.py [options]
Options:
  -h|--help:  show this help message
  -o|--output:  path of output files
  -l|--length:  a reference length file
  -f|--hicfile: path of high HiC matrix files (The matrix are named chr*.txt and chr*_chr*.txt)
  -a|--intra_binsize: the binsize of intra chromosome, default 50000
  -e|--inter_binsize: the binsize of inter chromosome, default 50000
  -w|--window:  local region window, default 10
  -c|--cutoff:  select SV events, default 0.6 (recommended parameters: [0.4, 0.5, 0.6]; low sequencing depth: set to 0.4 and high sequencing depth: set to 0.6 )
```
## Running the default example
```
python HiSV.py -o ../test_output 
               -l ../ref.len 
               -f ../test_matrix/ 
               -b 50000
               -w 10
               -c 0.5
```
## Output of HiSV
In the HiSV ouput directory, you will find
1. ``../test_output/HiSV_chr*_chr*_score.txt``  intermediate files used to identify the SV breakpoints 
2. ``../test_output/HiSV_intra_SV_result.txt``  the final integrated intra-chromosomal SVs
3. ``../test_output/HiSV_inter_SV_result.txt``  the final integrated inter-chromosomal SVs

```
cat HiSV_inter_SV_result.txt

chrom1  start_pos end_pos chrom2  start_pos end_pos
chr3  169100000 171100000 chr10 77900000  79700000
```

## Annotation of SV type
We provide type annotations for each structural variant.
```
Usage: SV_type.py [options]
Options:
  -h|--help:  show this help message
  -f|--hicfile: path of HiC matrix files
  -b|--binsize: binsize of HiC matrix, default 50000
  -i|--input:  path of result file from HiSV
```

## Determination of SV breakpoints
We provide accurate breakpoint identification when the Hi-C sample is a high-resolution data. 
```
Usage: determin_breakpoint.py [options]
Options:
  -h|--help:  show this help message
  -o|--output:  path of output files
  -f|--hicfile: path of high resolution HiC matrix files
  -b|--binsize: binsize of HiC matrix, default 50000
  -i|--input:  path of result file from HiSV
  -p|--position:  the position chromosomal of SV event (inter/intra)
```
