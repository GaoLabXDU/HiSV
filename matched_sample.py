import numpy as np
import math
from optparse import OptionParser
from scipy.stats import binom


def read_len(len_file, binSize):
    # read chr len file
    with open(len_file) as f:
        chr_id = []
        chr_len = []
        chr_bin_num = []
        for line in f:
            line_data = line.strip().split('\t')
            chr_id.append(line_data[0].strip('chr'))
            chr_len.append(int(line_data[1]))
            chr_bin_num.append(int(int(line_data[1]) / binSize))

    chr_bin_num = np.array(chr_bin_num)
    chr_id = np.array(chr_id)
    return chr_bin_num, chr_id


def log_ratio(case, control):
    for m in range(np.shape(case)[0]):
        for n in range(np.shape(case)[1]):
            if control[m][n] != 0 and case[m][n] != 0 :
                case[m][n] = math.log(case[m][n]/control[m][n], 2)
            else:
                case[m][n] = 0
    return case


def downsample(matrix, alpha):
    num = len(matrix[0])
    for i in range(num):
        for j in range(i, num):
            if matrix[i][j] != 0:
                cur_count = int(matrix[i][j])
                matrix[i][j] = binom.rvs(cur_count, alpha, size=1)
    return matrix


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c', '--case', default='/media/li/Data/GM12878/K562_remove_gap/', help='path of cancer matrix file (chr*.txt)')
    parser.add_option('-l', '--length', default='/home/li/SV/HiSV-main/test_data/ref.len', help='length file of reference')
    parser.add_option('-b', '--binsize', type='int', default=50000)
    parser.add_option('-n', '--normal',  default='/media/li/Data/GM12878/GM12878_remove_gap/', help="path of normal matrix file (chr*.txt)")
    parser.add_option('-o', '--output', default='/media/li/Data/GM12878/norm/', help='path of breakpoint file')

    (opts, args) = parser.parse_args()

    casefile = opts.case
    binSize = opts.binsize
    normalfile = opts.normal
    output = opts.output
    chr_len_file = opts.length

    chr_bin_num, chr_id = read_len(chr_len_file, binSize)

    for k in range(len(chr_id)):
        cur_chr = k
        num = chr_bin_num[cur_chr] + 1
        print(chr_id[cur_chr])

        case_file = casefile + "chr" + str(chr_id[cur_chr]) + ".txt"
        case_matrix = np.loadtxt(case_file)

        control_file = normalfile + "chr" + str(chr_id[cur_chr]) + ".txt"
        control_matrix = np.loadtxt(control_file)
        print(np.sum(case_matrix), np.sum(control_matrix))
        alpha = np.sum(case_matrix) / np.sum(control_matrix)
        print(alpha)
        if alpha > 1:
            case_matrix = downsample(case_matrix, alpha)
        else:
            control_matrix = downsample(control_matrix, alpha)

        norm_matrix = log_ratio(case_matrix, control_matrix)

        outfile = output + str(chr_id[cur_chr]) + ".txt"
        np.savetxt(outfile, norm_matrix, fmt="%i")
