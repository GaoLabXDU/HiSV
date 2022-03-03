import datetime
import numpy as np
import pysam
import pandas as pd
import math
import prox_tv as ptv
import os, sys
from optparse import OptionParser
from itertools import groupby
from collections import defaultdict


def read_reflen(chr_len_file):
    with open(chr_len_file) as f:
        chr_id = []
        chr_bin_num = []
        for line in f:
            line_data = line.strip().split('\t')
            chr_id.append(line_data[0])
            chr_bin_num.append(int(int(line_data[1]) / binSize))
    chr_bin_num = np.array(chr_bin_num)
    chr_start_pos = np.full(len(chr_bin_num)+1, 0)
    for i in range(1, len(chr_bin_num)+1):
        chr_start_pos[i] = sum(chr_bin_num[:i])
    return chr_id, chr_bin_num, chr_start_pos


def read_bam_file(filename, chrlist, matrix):
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


def local_saliency(matrix, win):
    d = int(win/2)
    N, M = len(matrix), len(matrix[0])
    nm_matrix = np.full((N, M), 0.0)
    Dist_matrix = np.full((win, win), 0.0)

    for m in range(win):
        for n in range(win):
            Dist_matrix[m][n] = np.sqrt(np.square(m) + np.square(n))

    it = np.nditer(matrix, flags=['multi_index'])
    while not it.finished:
        idx = it.multi_index
        data = matrix[idx]
        if data != 0:
            cur_matrix = matrix[idx[0]-d:idx[0]+d, idx[1]-d:idx[1]+d]
            cur_raw = cur_matrix.shape[0]
            cur_col = cur_matrix.shape[1]
            if cur_col == cur_raw == win and np.sum(cur_matrix) != 0:
                nm_matrix[idx] = 1 - math.exp(-np.mean(abs(data - cur_matrix) / (1 + Dist_matrix)))
            else:
                nm_matrix[idx] = 0
        it.iternext()
    return nm_matrix


def Calculating_diagonal_data(matrix):
    """
    # normalization matrix by diagonal to remove distance bias
    Calculating each diagonal mean and std
    """
    N, M = len(matrix), len(matrix[0])
    Diagonal_mean = np.full(M, 0.0)
    Diagonal_std = np.full(M, 0.0)
    std = []
    for d in range(N):
        intermediate = []
        c = d
        r = 0
        while r < N - d:
            intermediate.append(matrix[r][c])
            r += 1
            c += 1
        intermediate = np.array(intermediate)
        Diagonal_mean[d] = (np.mean(intermediate))
        Diagonal_std[d] = (np.std(intermediate))
    return Diagonal_mean, Diagonal_std


def Distance_normalization(matrix):

    """
    # normalization matrix by diagonal to remove distance bias
    norm_data = (data - mean_data) / mean_std
    """

    Diagonal_mean, Diagonal_std = Calculating_diagonal_data(matrix)
    N, M = len(matrix), len(matrix[0])
    for d in range(N):
        c = d
        r = 0
        while r < N - d:
            if Diagonal_std[d] == 0:
                matrix[r][c] = 0
            else:
                if matrix[r][c] - Diagonal_mean[d] < 0:
                    matrix[r][c] = 0
                else:
                    matrix[r][c] = (matrix[r][c] - Diagonal_mean[d]) / Diagonal_std[d]
            r += 1
            c += 1
    return matrix


def combine(lst):
    pos = (j - i for i, j in enumerate(lst))
    t = 0
    for i, els in groupby(pos):
        l = len(list(els))
        el = lst[t]
        t += l
        yield range(el, el+l)


def group(chrlist):
    last = chrlist[0]
    start = end = 0
    for n in chrlist[1:]:
        if n == last: # Part of the group, bump the end
            last = n
            end += 1
        else: # Not part of the group, yield current group and start a new
            yield start, end
            last = n
            start = end = end + 1
    yield start, end


def group_position(pos1list, pos2list):
    result = defaultdict(list)
    l1 = sorted(set(pos1list), key=pos1list.index)
    pos1group = list(combine(l1))
    for i in range(len(pos1group)):
        start_pos = pos1list.index(pos1group[i][0])
        end_pos = len(pos1list) - 1 - pos1list[::-1].index(pos1group[i][-1])
        # end_pos = pos1list.index(pos1group[i][-1])
        cur_list = pos2list[start_pos:end_pos+1]
        cur_list = sorted(cur_list)
        l2 = sorted(set(cur_list), key=cur_list.index)
        pos2group = list(combine(l2))
        if len(pos2group) > 1:
            for j in range(len(pos2group)):
                new_pos1 = []
                for m in range(len(pos1list)):
                    if pos2list[m] in pos2group[j]:
                        new_pos1.append(pos1list[m])
                new_pos1 = sorted(new_pos1)
                result['pos1_start'].append(new_pos1[0] * binSize)
                result['pos1_end'].append(new_pos1[-1] * binSize + binSize)
                result['pos2_start'].append(pos2group[j][0] * binSize)
                result['pos2_end'].append(pos2group[j][-1] * binSize + binSize)

        else:
            # print(pos1group[i][0], pos1group[i][-1], pos2group[0][0], pos2group[0][-1])
            result['pos1_start'].append(pos1group[i][0]*binSize)
            result['pos1_end'].append(pos1group[i][-1]*binSize + binSize)
            result['pos2_start'].append(pos2group[0][0]*binSize)
            result['pos2_end'].append(pos2group[0][-1]*binSize + binSize)

    return result

def group_1(data):
    last = data[0]
    start = end = 0
    for n in data[1:]:
        if n - last == 1:  # Part of the group, bump the end
            last = n
            end += 1
        else:  # Not part of the group, yield current group and start a new
            yield range(data[start], data[end] + 1)
            last = n
            start = end = end + 1
    # yield start, end
    yield range(data[start], data[end] + 1)


def group_inter_position(pos1list, pos2list):
    result = defaultdict(list)
    l1 = sorted(set(pos1list), key=pos1list.index)
    pos1group = list(group_1(l1))
    for i in range(len(pos1group)):
        start_pos = pos1list.index(pos1group[i][0])
        end_pos = len(pos1list) - 1 - pos1list[::-1].index(pos1group[i][-1])
        cur_pos2_list = pos2list[start_pos:end_pos + 1]
        cur_pos1_list = pos1list[start_pos:end_pos + 1]
        cur_pos2_list_sort = sorted(cur_pos2_list)
        l2 = sorted(set(cur_pos2_list_sort), key=cur_pos2_list_sort.index)
        pos2group = list(group_1(l2))
        if len(pos2group) > 1:
            for j in range(len(pos2group)):
                new_pos1 = []
                for k in range(len(pos2group[j])):
                    if pos2group[j][k] in cur_pos2_list:
                        pos = cur_pos2_list.index(pos2group[j][k])
                        new_pos1.append(cur_pos1_list[pos])
                '''
                for m in range(len(pos1list)):
                    if pos2list[m] in pos2group[j]:
                        new_pos1.append(pos1list[m])

                for m in range(len(cur_pos1_list)):
                    if pos2list[m] in cur_pos2_list:
                        new_pos1.append(cur_pos1_list[m])
                '''
                new_pos1 = sorted(new_pos1)
                result['pos1_start'].append(new_pos1[0] * binSize)
                result['pos1_end'].append(new_pos1[-1] * binSize + binSize)
                result['pos2_start'].append(pos2group[j][0] * binSize)
                result['pos2_end'].append(pos2group[j][-1] * binSize + binSize)

        else:
            result['pos1_start'].append(pos1group[i][0] * binSize)
            result['pos1_end'].append(pos1group[i][-1] * binSize + binSize)
            result['pos2_start'].append(pos2group[0][0] * binSize)
            result['pos2_end'].append(pos2group[0][-1] * binSize + binSize)
    return result


if __name__ == '__main__':

    starttime = datetime.datetime.now()
    parser = OptionParser()
    parser.add_option('-o', '--output', default='./output/', help='path of outfile')
    parser.add_option('-l', '--length', default='./ref.len', help='length file of reference')
    parser.add_option('-f', '--hicfile', default='./test.bam', help='hic data')
    parser.add_option('-b', '--binsize', type ='int', default=50000)
    parser.add_option('-w', '--window', type ='int', default=10)
    parser.add_option('-c', '--cutoff', type ='float', default=0.5)
    (opts, args) = parser.parse_args()

    output = opts.output
    if not os.path.exists(output):
        os.makedirs(output)
    
    chr_len_file = opts.length
    hicfile = opts.hicfile
    binSize = opts.binsize
    win = opts.window
    cutoff = opts.cutoff
    chr_id, chr_bin_num, chr_start_pos = read_reflen(chr_len_file)
    print(chr_bin_num)

    chr_count = len(chr_id)
    intra_result_file = output + 'intra_SV_result.txt'
    inter_result_file = output + 'inter_SV_result.txt'
    intra_chr = []
    intra_pos1 = []
    intra_pos2 = []
    inter_chr1 = []
    inter_chr2 = []
    inter_pos1 = []
    inter_pos2 = []
    for i in range(chr_count):
        for j in range(i, chr_count):
            if i == j:
                result_chr = pos1 = pos2 = []
                print('read Hi-C matrix: ', chr_id[i], chr_id[j])
                raw_num = chr_bin_num[i] + 1
                col_num = chr_bin_num[j] + 1
                chrlist = [chr_id[i], chr_id[j]]
                init_matrix = np.full((raw_num, col_num), 0.0)
                if 'matrix' in hicfile:
                    hic_matrix = np.loadtxt(hicfile)
                    contact_matrix = hic_matrix[chr_start_pos[i]:chr_start_pos[i+1], chr_start_pos[i]:chr_start_pos[i+1]]
                else:
                    contact_matrix = read_bam_file(bamfile, chrlist, init_matrix)
                il = np.tril_indices(raw_num)
                contact_matrix[il] = 0
                Dist_norm_matrix = Distance_normalization(contact_matrix)
                print('local saliency')
                local_sali_matrix = local_saliency(Dist_norm_matrix, win)
                tv_param = 0.2
                score = ptv.tv1_2d(local_sali_matrix, tv_param, n_threads=1, max_iters=0, method='dr')
                for m in range(raw_num):
                    for n in range(col_num):
                        if score[m][n] > cutoff:
                            intra_chr.append(chr_id[i])
                            intra_pos1.append(m)
                            intra_pos2.append(n)


            else:
                print(chr_id[i], chr_id[j])
                raw_num = chr_bin_num[i] + 1
                col_num = chr_bin_num[j] + 1
                chrlist = [chr_id[i], chr_id[j]]
                init_matrix = np.full((raw_num, col_num), 0.0)
                if 'matrix' in hicfile:
                    hic_matrix = np.loadtxt(hicfile)
                    contact_matrix = hic_matrix[chr_start_pos[i]:chr_start_pos[i+1], chr_start_pos[j]:chr_start_pos[j+1]]
                else:
                    contact_matrix = read_bam_file(bamfile, chrlist, init_matrix)

                print('local saliency')
                local_sali_matrix = local_saliency(contact_matrix, win)
                tv_param = 0.4
                score = ptv.tv1_2d(local_sali_matrix, tv_param, n_threads=1, max_iters=0, method='dr')
                for m in range(raw_num):
                    for n in range(col_num):
                        if score[m][n] > cutoff:
                            inter_chr1.append(chr_id[i])
                            inter_chr2.append(chr_id[j])
                            inter_pos1.append(m)
                            inter_pos2.append(n)

    # combine intra result
    if intra_chr:
        print(intra_chr)
        group_chr = list(group(intra_chr))
        for i in range(len(group_chr)):
            pos1list = intra_pos1[group_chr[i][0]:group_chr[i][1] + 1]
            pos2list = intra_pos2[group_chr[i][0]:group_chr[i][1] + 1]
            result = group_position(pos1list, pos2list)

            with open(intra_result_file, 'a') as out:
                for k in range(len(result['pos1_start'])):
                    out.write(str(intra_chr[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(result['pos1_end'][k])
                           + '\t' + str(result['pos2_start'][k]) + '\t' + str(result['pos2_end'][k]) + '\n')

    # combine inter result
    if inter_chr1:
        group_chr = list(group(inter_chr1))
        for i in range(len(group_chr)):
            chr2list = inter_chr2[group_chr[i][0]:group_chr[i][1]+1]
            cur_pos1 = inter_pos1[group_chr[i][0]:group_chr[i][1]+1]
            cur_pos2 = inter_pos2[group_chr[i][0]:group_chr[i][1]+1]
            group_chr2 = list(group(chr2list))
            for j in range(len(group_chr2)):
                pos1list = cur_pos1[group_chr2[j][0]:group_chr2[j][1] + 1]
                pos2list = cur_pos2[group_chr2[j][0]:group_chr2[j][1] + 1]
                result = group_inter_position(pos1list, pos2list)
                with open(inter_result_file, 'a') as out1:
                    for k in range(len(result['pos1_start'])):
                        out1.write(str(inter_chr1[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(result['pos1_end'][k])
                              + '\t' + str(chr2list[group_chr2[j][0]]) + '\t' + str(result['pos2_start'][k]) + '\t' + str(result['pos2_end'][k]) + '\n')


    endtime = datetime.datetime.now()
    print(endtime - starttime)

