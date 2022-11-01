import os
import numpy as np
import pysam


def read_reflen(chr_len_file, binSize):
    # read chrom length file
    # reture chr_id and chr_len/binsize(bin_num)
    with open(chr_len_file) as f:
        chr_id = []
        chr_bin_num = []
        for line in f:
            line_data = line.strip().split('\t')
            chr_id.append(line_data[0])
            chr_bin_num.append(int(int(line_data[1]) / binSize))

    chr_bin_num = np.array(chr_bin_num)
    return chr_id, chr_bin_num


def read_bam_file(filename, chrlist, matrix, binSize):
    bam = pysam.AlignmentFile(filename, "rb")
    for read in bam:
        if read.reference_name == chrlist[0] and read.next_reference_name == chrlist[1]:
            row = int(read.reference_start / binSize)
            col = int(read.next_reference_start / binSize)
            matrix[row][col] += 1
        elif read.reference_name == chrlist[1] and read.next_reference_name == chrlist[0]:
            row = int(read.next_reference_start / binSize)
            col = int(read.reference_start / binSize)
            matrix[row][col] += 1
    return matrix



def bamToMatrix(bamfile, ref_file, binSize, outdir, name):
    chr_id, chr_bin_num = read_reflen(ref_file, binSize)
    chr_count = len(chr_id)
    outdir_inter = os.path.join(outdir, "Inter_matrix")
    outdir_intra = os.path.join(outdir, "Intra_matrix")
    if not os.path.isdir(outdir_inter):
        os.mkdir(outdir_inter)
    if not os.path.isdir(outdir_intra):
        os.mkdir(outdir_intra)
    # i, j : chri and chrj
    for i in range(chr_count):
        for j in range(i, chr_count):
            if i == j:
                outputname = os.path.join(outdir_intra, name + '_%skb_%s_%s_matrix.txt' % (str(binSize//1000), chr_id[i], chr_id[j]))
            else:
                outputname = os.path.join(outdir_inter, name + '_%skb_%s_%s_matrix.txt' % (str(binSize//1000), chr_id[i], chr_id[j]))
            print(chr_id[i], chr_id[j])
            raw_num = chr_bin_num[i] + 1
            col_num = chr_bin_num[j] + 1
            chrlist = [chr_id[i], chr_id[j]]
            init_matrix = np.full((raw_num, col_num), 0)
            contact_matrix = read_bam_file(bamfile, chrlist, init_matrix, binSize)
            print(outputname)
            np.savetxt(outputname, contact_matrix, fmt='%.2f', delimiter='\t')
