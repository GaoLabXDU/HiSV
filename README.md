# HiSV
HiSV is a computational pipeline for structural variation detection from Hi-C data

Structural variations (SVs) play an essential role in the evolution of human genomes and are associated with cancer genetics and rare disease.

High-throughput chromosome capture (Hi-C) is a new technology that breaks through the limitations of the short-read-based whole-genome sequenc-ing (WGS) method and can be used to detect large-scale SVs. However, the identification of SVs from Hi-C data is still in the early stages with only a few methods available.

Here, we described HiSV (Hi-C for Structural Variation), a control-free computational pipeline-based salient object detection model to identify large-scale SVs from a Hi-C sample. 

# Required Packages
The following software must be installed on your machine:

python: tested with version 3.7

## Python dependencies:
* numpy (tested 1.21.6) 
* pysam (tested 0.21.0) (for hiccovert script)
* pandas (tested 1.3.5)
* pyBigWig (tested 0.3.22) (for SV_type script)
* prox_tv (tested 3.3.0) (for HiSV script)
* pytest (tested 7.4.0) 
* pygam (tested 0.8.0) (for SV_type script)
* scikit-learn (tested 1.0.2) (for determin_breakpoint script)
* cooler (tested 0.9.2) (for hiccovert script)
* hic-straw (tested 1.3.1) (for hiccovert script)
  
# Installation
HiSV and all the dependencies can be installed using pip:
```
conda create -n HiSV python=3.7
conda activate HiSV
pip install hisv -U
conda install hic-straw
```
# Overview
HiSV is distributed with four scripts. You can learn the basic usage of each script by typing command [-h] in a terminal window, where "command" is one of the following script names:
* __hiccovert__  
 Convert hic files in different formats to the matrix file format required by HiSV.
* __HiSV__  
 Calling SVs from Hi-C matrix data.
* __SV_type__  
 Type annotations for each structural variant.
* __determin_breakpoint__  
 Breakpoints for structural variants are more accurately identified when high-resolution Hi-C data are available.

# Tutorial
HiSV requires two input files, a tumor Hi-C data and a reference length file. HISV supports five types of Hi-C data as input including matrix file, bam file, .cool file, .hic file and fastq file. When Hi-C data is a bam file, .cool file, or hic file, you can use the script ``hiccovert`` to generate matrix file. When Hi-C data are fastq files, you can use the shell script ``fastqtobam.sh`` to generate bam file. This step requires the installation of BWA-MEM(https://github.com/lh3/bwa) and HiCExplorer(https://github.com/deeptools/HiCExplorer).  
HiSV results consist of two files including intra-chromosomal SV events and inter-chromosomal translocation events.   
Then, it provides ``SV_type`` and ``determin_breakpoint`` scripts to perform type annotations on SV events and more precise breakpoint identification respectively.
## Convert hic files in different formats to the matrix file format
Note: the input hic file format should be mcool, cool, hic or bam. The output Hi-C matrix file matches the number of chromosomes in the reference genome length file provided by the user
```
usage: hiccovert [-h] [--hic_file HIC_FILE] [--binsize BINSIZE]
                 [--format FORMAT] [--ref REF] [--cores CORES]
                 [--output OUTPUT] [--name NAME]

Convert Hi-C files in different formats into matrix form as HiSV input.

optional arguments:
  -h, --help           show this help message and exit
  --hic_file HIC_FILE  Hi-C file in different format. (default: None)
  --binsize BINSIZE    Bin size of final Hi-C contact matrix. (default: 50000)
  --format FORMAT      Format of input Hi-C file. This type needs to included
                       in [mcool, cool, hic, bam]. (default: None)
  --ref REF            Reference Genome length file (default: None)
  --cores CORES        The number of cores used for parallel computing
                       (default: None)
  --output OUTPUT      Output file path. (default: None)
  --name NAME          Sample name (default: None)
```
### Covert the Hi-C file in 'hic' format into matrix file
```
$ hiccovert --hic_file /mnt/d/HiSV_test_data/GSE63525_K562_combined.hic  --binsize 50000 --format hic  --ref /mnt/d/HiSV_test_data/ref.len  --output /mnt/d/HiSV_test_data/Matrix_data --name K562 --cores 20
```
### Covert the Hi-C file in 'bam' format into matrix file
```
$ hiccovert --hic_file /mnt/d/HiSV_test_data/K562_in_house_b38d5.nodup.bam  --binsize 50000 --format bam  --ref /mnt/d/HiSV_test_data/ref.len  --output /mnt/d/HiSV_test_data/Matrix_data --name K562 --cores 20
```
### Covert the Hi-C file in 'cool' format into matrix file
```
$ hiccovert --hic_file /mnt/d/HiSV_test_data/K562_50kb.cool  --binsize 50000 --format cool  --ref /mnt/d/HiSV_test_data/ref.len  --output /mnt/d/HiSV_test_data/Matrix_data --name K562 --cores 20
```
### Covert the Hi-C file in 'mcool' format into matrix file
```
$ hiccovert --hic_file /mnt/d/HiSV_test_data/K562.mcool  --binsize 50000 --format mcool  --ref /mnt/d/HiSV_test_data/ref.len  --output /mnt/d/HiSV_test_data/Matrix_data --name K562 --cores 20
```
### Example of reference Genome length file:
```
$ head ref.len
chr9    141213431
chr13   115169878
```

### The storage path of output matrix file is as follows:
```
--../HiSV_test_data/Matrix_data
---Intra_matrix
----K562_50Kb_chr9_chr9.matrix.txt
----K562_50Kb_chr13_chr13.matrix.txt
---Inter_matrix
----K562_50Kb_chr9_chr13.matrix.txt
```
## Calling SVs by HiSV
When the matrix files are prepared, the user can use the ``HiSV`` script to identify SVs from Hi-C matrix files. 
```
usage: HiSV [-h] [--inter_hic INTER_HIC] [--intra_hic INTRA_HIC]
            [--intra_binsize INTRA_BINSIZE] [--inter_binsize INTER_BINSIZE]
            [--window WINDOW] [--regularization REGULARIZATION]
            [--cutoff CUTOFF] [--ref REF] [--name NAME] [--cores CORES]
            [--output OUTPUT]

Identification of SVs based on Hi-C interaction matrix.

optional arguments:
  -h, --help            show this help message and exit
  --inter_hic INTER_HIC
                        inter chromosomal hi-c contact matrix. (default: None)
  --intra_hic INTRA_HIC
                        intra chromosomal hi-c contact matrix. (default: None)
  --intra_binsize INTRA_BINSIZE
                        Resolution of intra Hi-C contact matrix. (default:
                        50000)
  --inter_binsize INTER_BINSIZE
                        Resolution of inter Hi-C contact matrix. (default:
                        50000)
  --window WINDOW       Local region window. (default: 10)
  --regularization REGULARIZATION
                        The regularization parameter (default: 0.2)
  --cutoff CUTOFF       Threshold for filtering SV segments (default: 0.6)
  --ref REF             Reference Genome length file (default: None)
  --name NAME           Sample name (default: None)
  --cores CORES         The number of cores used for parallel computing
                        (default: None)
  --output OUTPUT       Output file path. (default: None)
```
### Note
By default, we set the local region window to 10, the regularization parameter to 0.2, and the cut-off values to 0.5 and 0.6 for Hi-C data with a sequencing depth of less than 100 million contacts and nearly 100 million contacts, respectively. We have found that these parameter settings achieve a good balance between sensitivity and precision for most test datasets.  
The local region is approcimately 1Mb could avoid the normal 3D organization having high salience. So we set the local region window to 10 for 50Kb resolution Hi-C data. When the resolution of Hi-C data is 100kb, we recommend setting window to 5.  
If your Hi-C data has a higher sequencing depth, we recommend increasing the cut-off value or the regularization parameter accordingly. 

Here, we provide an example Hi-C dataset to guide users through the whole process. The dataset contained two chromsomes, 9 and 13, and these two chromosomes in the K562 sample covered some of SVs. The command to run the test sample is：
```
$ HiSV --inter_hic /mnt/d/HiSV_test_data/Matrix_data/Inter_matrix --intra_hic /mnt/d/HiSV_test_data/Matrix_data/Intra_matrix --ref /mnt/d/HiSV_test_data/ref.len --name K562 --cores 20 --output /mnt/d/HiSV_test_data/result --cutoff 0.6
```

### Output of HiSV
In the HiSV ouput directory, you will find
1. ``../HiSV_test_data/result/K562_50Kb_chr*_chr*_HiSV_score.txt``  intermediate files used to identify the SV breakpoints 
2. ``../HiSV_test_data/result/HiSV_intra_SV_result.txt``  the integrated intra-chromosomal SVs
3. ``../HiSV_test_data/result/HiSV_inter_SV_result.txt``  the integrated inter-chromosomal SVs
The file "HiSV_intra_SV_result.txt" contains the intra-chromosomal SV events. Each column is named 'chrom1', 'start_pos', 'end_pos', 'start_pos' and 'end_pos'. Here, the breakpoint of each SV event is identified as a genome segment rather than a precise location.
```
$ head HiSV_intra_SV_result.txt
chr9    26550000        26800000        37400000        38450000
chr9    123550000       123600000       130900000       131150000
chr9    123550000       123600000       132400000       132550000
chr9    130900000       131100000       132400000       132550000
chr9    133700000       133900000       133800000       134000000
chr13   81050000        81300000        90400000        92550000
chr13   81050000        81150000        92950000        93350000
chr13   81050000        81100000        93850000        94050000
chr13   81050000        81250000        108500000       108750000
chr13   90450000        90500000        92950000        93350000
```
The file "HiSV_inter_SV_result.txt" contains the inter-chromosomal translocations. Each column is named 'chrom1', 'start_pos', 'end_pos', 'chrom2', 'start_pos' and 'end_pos'. Here, the breakpoint of each SV event is identified as a genome segment rather than a precise location.
```
$ head HiSV_inter_SV_result.txt
chr9    133550000       133750000       chr13   81000000        81550000
chr9    133500000       133650000       chr13   90400000        92500000
chr9    133600000       134050000       chr13   92850000        93350000
chr9    133550000       134000000       chr13   93850000        94200000
chr9    133550000       133950000       chr13   108350000       108950000
```
## Annotation of SV type
The next step is type annotations for each structural variant. Here, GC content, mappability and the number of restriction sites files are required as inputs to obtain a normalized coverage profile. These files are provided in /HiSV_test_data/normData, and some files with the extension of .gz need to be decompressed. Users can use the ``gzip -dv *`` command to decompress files.  
### Note  
This step only for intra-chromosomal SV events.
```
usage: SV_type [-h] [--hic_file HIC_FILE] [--gc_file GC_FILE]
               [--map_file MAP_FILE] [--RE_file RE_FILE] [--binsize BINSIZE]
               [--HiSV_result HISV_RESULT] [--name NAME]

Define the type of SV event recognized by each HiSV

optional arguments:
  -h, --help            show this help message and exit
  --hic_file HIC_FILE   Hi-C contact matrix file. (default: None)
  --gc_file GC_FILE     Hi-C contact matrix file. (default: None)
  --map_file MAP_FILE   Hi-C contact matrix file. (default: None)
  --RE_file RE_FILE     Hi-C contact matrix file. (default: None)
  --binsize BINSIZE     Bin size of final Hi-C contact matrix. (default:
                        50000)
  --HiSV_result HISV_RESULT
                        Result file from HiSV (default: None)
  --name NAME           Sample name (default: None)
```
The command to run the test sample is：
```
$ SV_type --hic_file /mnt/d/HiSV_test_data/Matrix_data/Intra_matrix --gc_file /mnt/d/HiSV_test_data/normData/hg19_1kb_GC.bw --map_file /mnt/d/HiSV_test_data/normData/hg19_mappability_100mer.1kb.bw --RE_file /mnt/d/HiSV_test_data/normData/hg19.HindIII.npz --binsize 50000 --HiSV_result /mnt/d/HiSV_test_data/result/HiSV_intra_SV_result.txt --name K562
```
### Output of SV_type
SV_type script outputs each SV type directly. For example:
```
________________________________________
SV_mean 105.09020572658713
all_mean 124.7138356396499
Num  1 deletion
________________________________________
SV_mean 147.46605838315728
all_mean 124.7138356396499
Num  2 duplication
______________________
```
SV_mean: the average coverage profile in SV region  
all_mean: the average coverage profile across the genome  
## Determination of SV breakpoints
When high-resolution Hi-C data is available, we provide the script to determine SV breakpoints more accurately.
```
usage: determin_breakpoint [-h] [--hic_file HIC_FILE] [--binsize BINSIZE]
                           [--HiSV_result HISV_RESULT] [--order_bin ORDER_BIN]
                           [--name NAME] [--output OUTPUT]

Determination of SV breakpoint

optional arguments:
  -h, --help            show this help message and exit
  --hic_file HIC_FILE   High resolution hi-C contact matrix file. (default:
                        None)
  --binsize BINSIZE     Bin size of high-resolution Hi-C contact matrix.
                        (default: 5000)
  --HiSV_result HISV_RESULT
                        Result file from HiSV (default: None)
  --order_bin ORDER_BIN
                        Bin size of Hi-C contact matrix for HiSV. (default:
                        50000)
  --name NAME           Sample name (default: None)
  --output OUTPUT       Output file path. (default: None)
```
Note: since the full storage of high-resolution Hi-C data occupies a lot of disk usage, sparse storage is used here to store Hi-C matrices, such as:
```
$ head K562_5kb_chr9_chr9_bed.txt
chrom_1 pos_1   chrom_2 pos_2   count
chr9    10000   chr9    10000   20
chr9    10000   chr9    15000   2
chr9    15000   chr9    15000   2
chr9    10000   chr9    20000   1
chr9    15000   chr9    20000   1
chr9    20000   chr9    20000   16
chr9    10000   chr9    25000   2
chr9    20000   chr9    25000   4
chr9    25000   chr9    25000   3
```
These high-reolution hi-c files are provided in ``/HiSV_test_data/Matrix_data/Intra_matrix`` and ``/HiSV_test_data/Matrix_data/Intra_matrix`` , and these files need to be decompressed. Users can use the ``gzip -dv *`` command to decompress files. The command to run the test sample is：
```
$ determin_breakpoint --hic_file /mnt/d/HiSV_test_data/Matrix_data/Intra_matrix --binsize 5000 --HiSV_result  /mnt/d/HiSV_test_data/result/HiSV_intra_SV_result.txt --name K562 --output /mnt/d/HiSV_test_data/result
```
### Output of determin_breakpoint
In the HiSV ouput directory, you will find
``../HiSV_test_data/result/HiSV_intra_SV_breakpoint.txt``
The file "HiSV_intra_SV_breakpoint.txt" contains the intra-chromosomal SV events. Each column is named 'chrom1', 'breakpoint1', 'chrom2' and 'breakpoint2'. Here, the breakpoint of each SV event is a precise location.
```
$ head HiSV_intra_SV_breakpoint.txt
chr9    26600000        chr9    38410000
chr9    123550000       chr9    130910000
chr9    123560000       chr9    132445000
chr9    131105000       chr9    132430000
chr9    133765000       chr9    133925000
chr13   81105000        chr13   90450000
chr13   81095000        chr13   93330000
chr13   81100000        chr13   93860000
chr13   81135000        chr13   108650000
chr13   90440000        chr13   92960000
```
