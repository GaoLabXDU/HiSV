#!/usr/bin/env python
import datetime
import numpy as np
import math
import prox_tv as ptv
import multiprocessing
import os, sys

import hisv.combineResult as combine_inter_result
import hisv.combineIntraResult as combine_intra_result
import argparse


def read_reflen(chr_len_file, inter_binSize, intra_binSize):
    """
    read chrom length file
    :param chr_len_file: reference genome length file
    :return: chromID
    """
    with open(chr_len_file) as f:
        chr_id = []
        chr_inter_count = []
        chr_intra_count = []
        for line in f:
            line_data = line.strip().split('\t')
            chr_id.append(line_data[0])
            chr_intra_count.append(int(line_data[1]) // intra_binSize + 1)
            chr_inter_count.append(int(line_data[1]) // inter_binSize + 1)
    return chr_id, chr_intra_count, chr_inter_count


def local_saliency(matrix, win):
    """
    # calculate local saliency map
    :param matrix: current Hi-C contact matrix
    :param win: 'Local region window.
    :return: saliency map
    """
    d = int(win / 2)
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
            cur_matrix = matrix[idx[0] - d:idx[0] + d, idx[1] - d:idx[1] + d]
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


def call_intra_sv(i, chr_id, chr_intra_count, chr_inter_count, intra_hicfile, inter_hicfile, win, reg, cutoff, name,
                  intra_binSize, inter_binSize, output):
    """
    Call SVs for each chromosome (including intra- and inter chromosomal)
    :param i: Chrom index
    :param chr_id: ChromID
    :param hicfile: current chromosomal hic matrix file
    :param win: local region window.
    :param cutoff: threshold for filtering SV segments
    """
    inter_chrom1 = []
    inter_chrom2 = []
    inter_pos1 = []
    inter_pos2 = []
    intra_chrom = []
    intra_pos1 = []
    intra_pos2 = []
    for j in range(i, len(chr_id)):
        chrlist = [chr_id[i], chr_id[j]]
        print("Read Hi-C matrix file: ", chrlist)

        if i == j:
            # intra_chromosomal
            score_result_file = os.path.join(output, name + '_%skb_%s_%s_HiSV_score.txt' % (str(intra_binSize // 1000),
                                                                                            chr_id[i], chr_id[j]))

            matrix_file = os.path.join(intra_hicfile, name + '_%skb_%s_%s_matrix.txt' % (str(intra_binSize // 1000),
                                                                                         chr_id[i], chr_id[j]))
            print(matrix_file)
            contact_matrix = np.loadtxt(matrix_file)
            if np.shape(contact_matrix)[0] != chr_intra_count[i]:
                print("The binsize setting does not match the matrix file.")
                sys.exit(-1)
            raw_num = np.shape(contact_matrix)[0]
            il = np.tril_indices(raw_num)
            contact_matrix[il] = 0
            # print('Z-score for distance')
            # z-score for distance normalization
            Dist_norm_matrix = Distance_normalization(contact_matrix)
            # print('local saliency')
            # calculating local saliency
            local_sali_matrix = local_saliency(Dist_norm_matrix, win)
            tv_param = reg
            score = ptv.tv1_2d(local_sali_matrix, tv_param, n_threads=1, max_iters=0, method='dr')
            np.savetxt(score_result_file, score, fmt="%1.2f")
            raw_num = col_num = np.shape(score)[0]

            for m in range(raw_num):
                for n in range(col_num):
                    if score[m][n] > cutoff:
                        intra_chrom.append(chr_id[i])
                        intra_pos1.append(m)
                        intra_pos2.append(n)

        else:
            # intrer_chromosomal
            score_result_file = os.path.join(output, name + '_%skb_%s_%s_HiSV_score.txt' % (str(inter_binSize // 1000),
                                                                                            chr_id[i], chr_id[j]))
            matrix_file = os.path.join(inter_hicfile, name + '_%skb_%s_%s_matrix.txt' % (str(inter_binSize // 1000),
                                                                                         chr_id[i], chr_id[j]))
            print(matrix_file)
            contact_matrix = np.loadtxt(matrix_file)
            if np.shape(contact_matrix)[0] != chr_inter_count[i] or np.shape(contact_matrix)[1] != chr_inter_count[j]:
                print("The binsize setting does not match the matrix file.")
                sys.exit(-1)
            raw_num = np.shape(contact_matrix)[0]
            col_num = np.shape(contact_matrix)[1]

            # calculating local saliency
            local_sali_matrix = local_saliency(contact_matrix, 10)
            tv_param = reg
            score = ptv.tv1_2d(local_sali_matrix, tv_param, n_threads=1, max_iters=0, method='dr')
            np.savetxt(score_result_file, score, fmt="%1.2f")
            for m in range(raw_num):
                for n in range(col_num):
                    if score[m][n] > cutoff:
                        inter_chrom1.append(chr_id[i])
                        inter_chrom2.append(chr_id[j])
                        inter_pos1.append(m)
                        inter_pos2.append(n)
    return intra_chrom, intra_pos1, intra_pos2, inter_chrom1, inter_chrom2, inter_pos1, inter_pos2


def getargs():
    ## Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(description='''Identification of SVs based on Hi-C interaction matrix.''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--inter_hic', help='inter chromosomal hi-c contact matrix.')
    parser.add_argument('--intra_hic', help='intra chromosomal hi-c contact matrix.')
    parser.add_argument('--intra_binsize', type=int, default=50000, help='Resolution of intra Hi-C contact matrix.')
    parser.add_argument('--inter_binsize', type=int, default=50000, help='Resolution of inter Hi-C contact matrix.')
    parser.add_argument('--window', type=int, default=10, help='''Local region window.''')
    parser.add_argument('--regularization', type=float, default=0.2, help='The regularization parameter')
    parser.add_argument('--cutoff', type=float, default=0.6, help='''Threshold for filtering SV segments''')
    parser.add_argument('--ref', help='Reference Genome length file')
    parser.add_argument('--name', help='Sample name')
    parser.add_argument('--cores', type=int, help='The number of cores used for parallel computing')
    parser.add_argument('--output', help='''Output file path.''')

    ## Parse the command-line arguments
    commands = sys.argv[1:]
    if not commands:
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands


def run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name, output):
    if not os.path.exists(output):
        os.makedirs(output)

    if cutoff > 1 or cutoff < 0:
        print("Cut-off parameter should be in the range 0-1.")
        # sys.exit(-1)
        return "Parameter Error."

    if reg > 1 or reg < 0:
        print("Regularization parameter should be in the range 0-1.")
        return "Parameter Error."

    if not isinstance(win, int) or win < 0:
        print("The Win parameter should be an integer greater than 0.")
        # sys.exit(-1)
        return "Parameter Error."

    if not os.path.exists(intra_hicfile):
        print("Intra_hicfile does not exist")
        # sys.exit(-1)
        return "File does not exist."

    if not os.path.exists(inter_hicfile):
        print("Inter_hicfile does not exist")
        # sys.exit(-1)
        return "File does not exist."

    if not os.path.exists(ref_file):
        print("reference file does not exist")
        # sys.exit(-1)
        return "File does not exist."

    intra_result_file = os.path.join(output, 'HiSV_intra_SV_result.txt')
    inter_result_file = os.path.join(output, 'HiSV_inter_SV_result.txt')
    inter_chrom1 = []
    inter_chrom2 = []
    inter_pos1 = []
    inter_pos2 = []
    intra_chrom = []
    intra_pos1 = []
    intra_pos2 = []

    # read chromosomal length file
    chr_id, chr_intra_count, chr_inter_count = read_reflen(ref_file, inter_binSize=inter_binSize,
                                                            intra_binSize=intra_binSize)

    # Determine whether the hi-c matrix file is consistent with the reference file
    for i in range(len(chr_id)):
        for j in range(i, len(chr_id)):
            if i == j:
                matrix_file = os.path.join(intra_hicfile, name + '_%skb_%s_%s_matrix.txt' % (str(intra_binSize // 1000),
                                                                                             chr_id[i], chr_id[j]))
            else:
                matrix_file = os.path.join(inter_hicfile, name + '_%skb_%s_%s_matrix.txt' % (str(inter_binSize // 1000),
                                                                                             chr_id[i], chr_id[j]))
            if not os.path.exists(matrix_file):
                print("The required matrix file should exist in the current folder.")
                return "Matrix file Error."

    # the number of chrom
    chr_count = len(chr_id)
    # multi-processing
    p = multiprocessing.Pool(processes=n_cores)
    # i, chr_id, intra_hicfile, inter_hicfile, win, cutoff, name, binSize
    async_results = [p.apply_async(call_intra_sv, args=(i, chr_id, chr_intra_count, chr_inter_count, intra_hicfile,
                                                        inter_hicfile, win, reg, cutoff, name, intra_binSize,
                                                        inter_binSize, output)) for i in range(chr_count)]

    p.close()
    p.join()

    for ar in async_results:
        cur_chrom_result = ar.get()
        # intra_chrom, intra_pos1, intra_pos2, inter_chrom1, inter_chrom2, inter_pos1, inter_pos2
        intra_chrom.extend(cur_chrom_result[0])
        intra_pos1.extend(cur_chrom_result[1])
        intra_pos2.extend(cur_chrom_result[2])
        inter_chrom1.extend(cur_chrom_result[3])
        inter_chrom2.extend(cur_chrom_result[4])
        inter_pos1.extend(cur_chrom_result[5])
        inter_pos2.extend(cur_chrom_result[6])

    # combine intra-chromosomal SV region
    if intra_pos1:
        group_chr = list(combine_intra_result.group(intra_chrom))
        for i in range(len(group_chr)):
            pos1list = intra_pos1[group_chr[i][0]:group_chr[i][1] + 1]
            pos2list = intra_pos2[group_chr[i][0]:group_chr[i][1] + 1]
            result = combine_intra_result.group_position(pos1list, pos2list, intra_binSize)
            with open(intra_result_file, 'a') as out1:
                for k in range(len(result['pos1_start'])):
                    out1.write(
                        str(intra_chrom[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(
                            result['pos1_end'][k])
                        + '\t' + str(result['pos2_start'][k]) + '\t' + str(result['pos2_end'][k]) + '\n')

    # combine inter-chromosomal SV region
    if inter_pos1:
        group_chr = list(combine_inter_result.group(inter_chrom1))
        for i in range(len(group_chr)):
            cur_chr2 = inter_chrom2[group_chr[i][0]:group_chr[i][1] + 1]
            no_order_pos1 = inter_pos1[group_chr[i][0]:group_chr[i][1] + 1]
            no_order_pos2 = inter_pos2[group_chr[i][0]:group_chr[i][1] + 1]
            sort_id = sorted(range(len(cur_chr2)), key=lambda k: cur_chr2[k])
            chr2list = [cur_chr2[i] for i in sort_id]
            cur_pos1 = [no_order_pos1[i] for i in sort_id]
            cur_pos2 = [no_order_pos2[i] for i in sort_id]
            group_chr2 = list(combine_inter_result.group(chr2list))

            for j in range(len(group_chr2)):
                pos1list = cur_pos1[group_chr2[j][0]:group_chr2[j][1] + 1]
                pos2list = cur_pos2[group_chr2[j][0]:group_chr2[j][1] + 1]
                result = combine_inter_result.group_position(pos1list, pos2list, inter_binSize)
                with open(inter_result_file, 'a') as out1:
                    for k in range(len(result['pos1_start'])):
                        out1.write(
                            str(inter_chrom1[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(
                                result['pos1_end'][k])
                            + '\t' + str(chr2list[group_chr2[j][0]]) + '\t' + str(result['pos2_start'][k]) + '\t' + str(
                                result['pos2_end'][k]) + '\n')

    return "Finish."


if __name__ == '__main__':

    start_t = datetime.datetime.now()
    args, commands = getargs()

    output = args.output
    if not os.path.exists(output):
        os.makedirs(output)
    ref_file = args.ref
    inter_binSize = args.inter_binsize
    intra_binSize = args.intra_binsize
    win = args.window
    cutoff = args.cutoff
    reg = args.regularization
    n_cores = args.cores
    intra_hicfile = args.intra_hic
    inter_hicfile = args.inter_hic
    name = args.name
    # ref_file, inter_binSize, intra_binSize, win, cutoff, n_cores, intra_hicfile, inter_hicfile, name, output
    run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name, output)
    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("running time: " + "{:.2f}".format(elapsed_sec) + " seconds")
