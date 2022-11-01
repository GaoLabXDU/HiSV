import numpy as np
import os
from collections import defaultdict
import seaborn
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from sklearn.decomposition import PCA
from optparse import OptionParser

from itertools import groupby

"""
Recalculate SV breakpoint in high resolution matrix
Reference
https://github.com/XiaoTaoWang/NeoLoopFinder
"""


def get_SV_intra_info(SV_file):
    SV_line = open(SV_file)
    SV_info = defaultdict(list)
    for line in SV_line:
        line = line.strip().split('\t')
        # if line[0] in chrlist:
        SV_info['chrom1'].append(line[0])
        SV_info['bp1_start'].append(int(line[1])-50000)
        SV_info['bp1_end'].append(int(line[2])+50000)
        SV_info['chrom2'].append(line[0])
        SV_info['bp2_start'].append(int(line[3])-50000)
        SV_info['bp2_end'].append(int(line[4])+50000)
    return SV_info


def get_SV_inter_info(SV_file):
    SV_line = open(SV_file)
    SV_info = defaultdict(list)
    for line in SV_line:
        line = line.strip().split('\t')
        # if line[0] in chrlist:
        SV_info['chrom1'].append(line[0])
        SV_info['bp1_start'].append(int(line[1])-50000)
        SV_info['bp1_end'].append(int(line[2])+50000)
        SV_info['chrom2'].append(line[3])
        SV_info['bp2_start'].append(int(line[4])-50000)
        SV_info['bp2_end'].append(int(line[5])+50000)
    return SV_info


def group(data):
    pos = []
    for i in range(len(data)-1):
        if data[i] != data[i+1]:
            pos.append(i)
    return pos


def recalculate_breakpoint(matrix):
    # row
    differ_value = []
    corr = gaussian_filter(np.corrcoef(matrix, rowvar=True), sigma=1)

    corr[np.isnan(corr)] = 0
    pca = PCA(n_components=3, whiten=True)
    pca1 = pca.fit_transform(corr)[:, 0]
    # loci_1 = np.where(pca1 >= 0)[0]

    loci_1 = np.sign(pca1)
    poslist = group(loci_1)
    if poslist:
        for cur_pos in poslist:
            differ_value.append(abs(pca1[cur_pos]-pca1[cur_pos+1]))
        index = differ_value.index(max(differ_value))
        up_row = poslist[index]

    else:
        up_row = 0

    # column
    differ_value = []
    corr = gaussian_filter(np.corrcoef(matrix, rowvar=False), sigma=1)
    corr[np.isnan(corr)] = 0
    pca = PCA(n_components=3, whiten=True)
    pca2 = pca.fit_transform(corr)[:, 0]

    loci_2 = np.sign(pca2)
    poslist = group(loci_2)
    if poslist:
        for cur_pos in poslist:
            differ_value.append(abs(pca2[cur_pos]-pca2[cur_pos+1]))
        index = differ_value.index(max(differ_value))
        up_col = poslist[index]
    else:
        up_col = 0
    return up_row, up_col


def get_SV_matrix(SV_region, matrix_file, binSize):
    row_num = (SV_region['bp1_end'] - SV_region['bp1_start']) // binSize + 1
    col_num = (SV_region['bp2_end'] - SV_region['bp2_start']) // binSize + 1
    Init_matrix = np.full((row_num, col_num), 0)
    bed_list = open(matrix_file)
    bed_list.readline()
    for line in bed_list:
        line = line.strip().split(' ')
        pos1 = int(line[1])
        pos2 = int(line[2])
        if SV_region['bp1_start'] <= pos1 <= SV_region['bp1_end'] and SV_region['bp2_start'] <= pos2 <= SV_region['bp2_end']:

            cur_pos1 = (pos1 - SV_region['bp1_start']) // binSize
            cur_pos2 = (pos2 - SV_region['bp2_start']) // binSize
            Init_matrix[cur_pos1][cur_pos2] = int(line[3])
    return Init_matrix



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-f', '--hicfile', default='/media/li/Data/K562/5kb/intra/', help='path of hic matrix file')
    parser.add_option('-b', '--binsize', type='int', default=5000)
    parser.add_option('-i', '--input',  default='/media/li/Data/GM12878/K562_matrix/K562_combine_result.txt', help="path of result file from HiSV")
    parser.add_option('-p', '--position',  default='intra', help="intra/inter")
    parser.add_option('-o', '--output', default='/media/li/Data/GM12878/K562_matrix', help='path of breakpoint file')

    (opts, args) = parser.parse_args()

    hicfile = opts.hicfile
    binSize = opts.binsize
    SV_file = opts.input
    position = opts.position
    output = opts.output
    output_file = output + "/HiSV_" + str(position) + "_SV_breakpoint_result.txt"
    if position == "inter":
        SV_info = get_SV_inter_info(SV_file)
    else:
        SV_info = get_SV_intra_info(SV_file)

    with open(output_file, 'a') as out:
        for i in range(1):
            chrom1 = SV_info['chrom1'][i]
            chrom2 = SV_info['chrom2'][i]
            matrix_file = hicfile + str(chrom1) + '_' + str(chrom2) + '.txt'
            print(matrix_file)
            SV_ragion = {}
            SV_ragion['bp1_start'] = SV_info['bp1_start'][i]
            SV_ragion['bp1_end'] = SV_info['bp1_end'][i]
            SV_ragion['bp2_start'] = SV_info['bp2_start'][i]
            SV_ragion['bp2_end'] = SV_info['bp2_end'][i]
            matrix = get_SV_matrix(SV_ragion, matrix_file, binSize)
            row_pos, column_pos = recalculate_breakpoint(matrix)
            start_pos = row_pos * binSize + SV_info['bp1_start'][i]
            end_pos = column_pos * binSize + SV_info['bp2_start'][i]
            print(row_pos, column_pos, start_pos, end_pos)
            out.write(str(chrom1) + '\t' + str(start_pos) + '\t' + str(chrom2) + '\t' + str(end_pos) + '\n')

