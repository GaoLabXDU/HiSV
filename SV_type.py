from collections import defaultdict
import numpy as np
from iced import normalization
from optparse import OptionParser


def read_file(resultfile):
    result = defaultdict(list)
    with open(resultfile) as f:
        for line in f:
            line_data = line.strip().split('\t')
            result['chr'].append(line_data[0])
            result['pos1_start'].append(int(line_data[1]))
            result['pos1_end'].append(int(line_data[2]))
            result['pos2_start'].append(int(line_data[3]))
            result['pos2_end'].append(int(line_data[4]))
    return result


def call_slope(data):
    order = 1
    index = [i for i in range(1, len(data) + 1)]
    coeffs = np.polyfit(index, list(data), order)
    slope = coeffs[-2]
    return slope


def cal_slope(matrix):
    data_1 = matrix[:, 0]
    data_2 = matrix[:, -1]
    slope1 = call_slope(data_1)
    slope2 = call_slope(data_2)
    return slope1, slope2


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-f', '--hicfile', default='/media/li/Data/GM12878/K562_remove_gap/', help='path of hic matrix file (chr*.txt)')
    parser.add_option('-b', '--binsize', type='int', default=50000)
    parser.add_option('-i', '--input',  default='/media/li/Data/GM12878/K562_matrix/K562_combine_result.txt', help="path of result file from HiSV")
    (opts, args) = parser.parse_args()


    hicfile = opts.hicfile
    binSize = opts.binsize
    input_file = opts.input

    # read result file from HiSV
    result = read_file(input_file)
    for i in range(len(result['chr'])):
        print(result['chr'][i])
        # read matrix file
        matrix_file = hicfile + result['chr'][i] + ".txt"
        contact_matrix = np.loadtxt(matrix_file)
        # normalization
        normed_matrix = normalization.ICE_normalization(contact_matrix)
        raw_num = len(normed_matrix[0])
        count_mean = np.full(raw_num, 0.0)
        for k in range(raw_num):
            count_mean[k] = np.mean(normed_matrix[k])
        # average coverage profile
        all_mean = np.mean(count_mean)
        print(all_mean)

        start = int(result['pos1_start'][i] / binSize)
        end = int(result['pos2_end'][i] / binSize)
        print('SV_mean', np.mean(count_mean[start:end]))
        print('all_mean', np.mean(count_mean))
        if np.mean(count_mean[start:end]) > all_mean * 1.1:
            print("Num ", str(i+1), 'duplication')
        elif np.mean(count_mean[start:end]) < all_mean / 1.1:
            print("Num ", str(i+1), 'deletion')
        else:
            start_1 = int(result['pos1_start'][i] / binSize)
            end_1 = int(result['pos1_end'][i] / binSize)
            start_2 = int(result['pos2_start'][i] / binSize)
            end_2 = int(result['pos2_end'][i] / binSize)
            cur_matrix = normed_matrix[start_1:end_1, start_2:end_2]
            slope1, slope2 = cal_slope(cur_matrix)
            print(slope1, slope2)
            if slope1 >= 0 and slope2 >= 0:
                print("Num ", str(i+1), 'unbalanced_translocation')
            elif slope1 <= 0 and slope2 <= 0:
                print("Num ", str(i+1), 'unbalanced_translocation')
            elif slope1 < 0 and slope2 > 0:
                print("Num ", str(i+1), 'inversion')
            else:
                print("Num ", str(i+1), 'balanced_translocation')


