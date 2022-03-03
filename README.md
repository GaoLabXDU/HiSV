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
HiSV requires two input files, a tumor Hi-C data and a reference length file. HISV supports three types of Hi-C data as input including matrix file, bam file and fastq file. When Hi-C data are fastq files, you can use the shell script fastqtobam.sh to generate bam file. This step requires the installation of BWA-MEM(https://github.com/lh3/bwa) and HiCExplorer(https://github.com/deeptools/HiCExplorer). HiSV results consist of two files including intra-chromosomal SV events and inter-chromosomal translocation events. 

## Ruuning command
```
Usage: HiSV.py [options]
Options:
  -h|--help:  show this help message
  -o|--output:  path of output files
  -l|--length:  a reference length file
  -f|--hicfile: Hi-C data (_.matrix or _.bam)
  -b|--binsize: binsize of HiC matrix, default 50000
  -w|--window:  local region window, default 10
  -c|--cutoff:  select SV events, default 0.5
```
## Running the default example
```
python HiSV.py -o ../test_output 
               -l ../ref.len 
               -f ../test.bam 
               -b 50000
               -w 10
               -c 0.5
```

 




