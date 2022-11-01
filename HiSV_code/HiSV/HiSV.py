import datetime
import numpy as np
import math
import prox_tv as ptv
import os, sys
from optparse import OptionParser
import combine_result as combine_inter_result
import combine_intra_result
import determin_breakpoint


def read_reflen(chr_len_file):
    with open(chr_len_file) as f:
        chr_id = []
        chr_bin_num = []
        for line in f:
            line_data = line.strip().split('\t')
            chr_id.append(line_data[0])
            chr_bin_num.append(int(int(line_data[1]) / binSize))

    chr_bin_num = np.array(chr_bin_num)
    return chr_id, chr_bin_num


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


if __name__ == '__main__':

    starttime = datetime.datetime.now()
    parser = OptionParser()
    parser.add_option('-o', '--output', default='./output', help='path of out file')
    parser.add_option('-l', '--length', default='./ref.len', help='length file of reference')
    parser.add_option('-f', '--hicfile', default='./test.bam', help='path of hic matrix file')
    parser.add_option('-b', '--binsize', type='int', default=50000)
    parser.add_option('-w', '--window', type='int', default=10)
    parser.add_option('-c', '--cutoff', type='float', default=0.3)
    (opts, args) = parser.parse_args()

    output = opts.output
    if not os.path.exists(output):
        os.makedirs(output)

    chr_len_file = opts.length
    hicfile = opts.hicfile
    binSize = opts.binsize
    win = opts.window
    cutoff = opts.cutoff
    intra_result_file = output + '/HiSV_intra_SV_result.txt'
    inter_result_file = output + '/HiSV_inter_SV_result.txt'
    inter_chrom1 = []
    inter_chrom2 = []
    inter_pos1 = []
    inter_pos2 = []
    intra_chrom = []
    intra_pos1 = []
    intra_pos2 = []

    # read chromsomal length file
    chr_id, chr_bin_num = read_reflen(chr_len_file)

    chr_count = len(chr_id)

    print("Read Hi-C matrix file...")
    print(chr_id, chr_bin_num)
    for i in range(chr_count):
        for j in range(i, chr_count):
            chrlist = [chr_id[i], chr_id[j]]
            print("Read Hi-C matrix file: ", chrlist)
            score_result_file = output + '/HiSV_' + chr_id[i] + '_' + chr_id[j] + "_score.txt"
            if i == j:
                # intra_chromosomal
                matrix_file = hicfile + str(chr_id[i]) + ".txt"
                print(matrix_file)
                contact_matrix = np.loadtxt(matrix_file)
                raw_num = col_num = np.shape(contact_matrix)[0]

                il = np.tril_indices(raw_num)
                contact_matrix[il] = 0

                print('Z-score for distance')
                Dist_norm_matrix = Distance_normalization(contact_matrix)

                print('local saliency')
                local_sali_matrix = local_saliency(Dist_norm_matrix, 10)

                tv_param = 0.2
                score = ptv.tv1_2d(local_sali_matrix, tv_param, n_threads=1, max_iters=0, method='dr')

            # with open(outfile, 'a') as out:
                # outfile = "/media/li/Data/GM12878/K562_remove_gap/test_result/" + "HiSV_score_matrix_" + str(chr_id[i]) + '_' + str(chr_id[j]) + ".txt"
                np.savetxt(score_result_file, score, fmt="%1.2f")
                for m in range(raw_num):
                    for n in range(col_num):
                        if score[m][n] > cutoff:
                            intra_chrom.append(chr_id[i])
                            intra_pos1.append(m)
                            intra_pos2.append(n)
                            # out.write(str(chr_id[i]) + '\t' + str(m*binSize) + '\t' + str(n*binSize) + '\n')
                            # out.write(str(chr_id[i]) + '\t' + str(m*binSize) + '\t' + str(chr_id[j]) + '\t' + str(n*binSize) + '\n')
                            # print(i + 1, m * binsize, j + 1, n * binsize, score_matrix[m][n])
            else:
                # intrer_chromosomal
                matrix_file = hicfile + str(chr_id[i]) + "_" + str(chr_id[j]) +  ".txt"
                print(matrix_file)
                contact_matrix = np.loadtxt(matrix_file)
                raw_num = np.shape(contact_matrix)[0]
                col_num = np.shape(contact_matrix)[1]

                print('local saliency')
                local_sali_matrix = local_saliency(contact_matrix, 10)
                tv_param = 0.2
                score = ptv.tv1_2d(local_sali_matrix, tv_param, n_threads=1, max_iters=0, method='dr')
                np.savetxt(score_result_file, score, fmt="%1.2f")
                for m in range(raw_num):
                    for n in range(col_num):
                        if score[m][n] > cutoff:
                            inter_chrom1.append(chr_id[i])
                            inter_chrom2.append(chr_id[j])
                            inter_pos1.append(m)
                            inter_pos2.append(n)

    # combine intra-chromosomal SV region
    group_chr = list(combine_intra_result.group(intra_chrom))
    for i in range(len(group_chr)):
        pos1list = intra_pos1[group_chr[i][0]:group_chr[i][1] + 1]
        pos2list = intra_pos2[group_chr[i][0]:group_chr[i][1] + 1]
        result = combine_intra_result.group_position(pos1list, pos2list, binSize)
        with open(intra_result_file, 'a') as out1:
            for k in range(len(result['pos1_start'])):
                out1.write(
                    str(intra_chrom[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(result['pos1_end'][k])
                    + '\t' + str(result['pos2_start'][k]) + '\t' + str(result['pos2_end'][k]) + '\n')

    # combine inter-chromosomal SV region
    group_chr = list(combine_inter_result.group(inter_chrom1))
    for i in range(len(group_chr)):
        cur_chr2 = inter_chrom2[group_chr[i][0]:group_chr[i][1]+1]
        no_order_pos1 = inter_pos1[group_chr[i][0]:group_chr[i][1]+1]
        no_order_pos2 = inter_pos2[group_chr[i][0]:group_chr[i][1]+1]
        sort_id = sorted(range(len(cur_chr2)), key=lambda k: cur_chr2[k])
        chr2list = [cur_chr2[i] for i in sort_id]
        cur_pos1 = [no_order_pos1[i] for i in sort_id]
        cur_pos2 = [no_order_pos2[i] for i in sort_id]
        group_chr2 = list(combine_inter_result.group(chr2list))

        for j in range(len(group_chr2)):
            pos1list = cur_pos1[group_chr2[j][0]:group_chr2[j][1] + 1]
            pos2list = cur_pos2[group_chr2[j][0]:group_chr2[j][1] + 1]
            result = combine_inter_result.group_position(pos1list, pos2list, binSize)
            with open(inter_result_file, 'a') as out1:
                for k in range(len(result['pos1_start'])):
                    out1.write(str(inter_chrom1[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(result['pos1_end'][k])
                              + '\t' + str(chr2list[group_chr2[j][0]]) + '\t' + str(result['pos2_start'][k]) + '\t' + str(result['pos2_end'][k]) + '\n')







