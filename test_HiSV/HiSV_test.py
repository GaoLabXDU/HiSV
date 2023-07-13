import pytest
import pandas as pd
from HiSV import run_HiSV
import os


def read_intra_groundtruth(filename):
    intra_GT = pd.DataFrame(columns=['chrom', 'start', 'end'])
    with open(filename) as f:
        for line in f:
            line_data = line.strip().split('\t')
            intra_GT = intra_GT.append(
                {'chrom': line_data[0], 'start': int(line_data[1]), 'end': int(line_data[2])}, ignore_index=True)
    return intra_GT


def sta_intra_result(intra_GT, intra_filename):
    count = 0
    result_count = 0
    with open(intra_filename) as f:
        for line in f:
            line_data = line.strip().split('\t')
            result_count += 1
            result_chr = line_data[0]
            result_start = int(line_data[1])
            result_end = int(line_data[4])
            for i in intra_GT.index:
                GT_row = intra_GT.loc[i]
                # 'chrom', 'start', 'end'
                GT_chr = GT_row['chrom']
                GT_start = GT_row['start']
                GT_end = GT_row['end']

                if result_chr == GT_chr:
                    if GT_start <= result_start <= GT_end <= result_end:
                        if (GT_end - result_start) >= (result_end - result_start) * 0.75:
                            count += 1
                            break
                    elif GT_end >= result_end >= GT_start >= result_start:
                        if (result_end - GT_start) >= (result_end - result_start) * 0.75:
                            count += 1
                            break
                    elif result_start <= GT_start <= result_end and result_start <= GT_end <= result_end:
                        if (GT_end - GT_start) >= (result_end - result_start) * 0.75:
                            count += 1
                            break

    f1_score = 2 * (count / len(intra_GT)) * (count / result_count) / (count / len(intra_GT) + (count / result_count))
    return f1_score


def read_inter_groundtruth(filename):
    inter_GT = pd.DataFrame(columns=['chrom1', 'pos1', 'chrom2', 'pos2'])
    with open(filename) as f:
        for line in f:
            line_data = line.strip().split('\t')
            inter_GT = inter_GT.append(
                {'chrom1': line_data[0], 'pos1': int(line_data[1]),
                 'chrom2': line_data[2], 'pos2': int(line_data[3])}, ignore_index=True)
    return inter_GT


def sta_inter_result(inter_GT, inter_filename, binSize):
    count = 0
    result_count = 0
    with open(inter_filename) as f:
        for line in f:
            line_data = line.strip().split('\t')
            result_count += 1
            result_chr_1 = line_data[0]
            result_start_1 = int(line_data[1]) - binSize
            result_end_1 = int(line_data[2]) + binSize
            result_chr_2 = line_data[3]
            result_start_2 = int(line_data[4]) - binSize
            result_end_2 = int(line_data[5]) + binSize
            for i in inter_GT.index:
                GT_row = inter_GT.loc[i]
                # 'chrom1', 'pos1', 'chrom2', 'pos2'
                GT_chr1 = GT_row['chrom1']
                GT_chr2 = GT_row['chrom2']
                GT_pos1 = GT_row['pos1']
                GT_pos2 = GT_row['pos2']
                if result_chr_1 == GT_chr1 and result_chr_2 == GT_chr2:
                    if result_start_1 <= GT_pos1 <= result_end_1 and result_start_2 <= GT_pos2 <= \
                            result_end_2:
                        count += 1

                elif result_chr_1 == GT_chr2 and result_chr_2 == GT_chr1:
                    if result_start_2 <= GT_pos1 <= result_end_2 and result_start_1 <= GT_pos2 <= \
                            result_end_1:
                        count += 1
    f1_score = 2 * (count / len(inter_GT)) * (count / result_count) / (count / len(inter_GT) + (count / result_count))
    return f1_score


def test_fileExist():
    # file not exist
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "ref_1.len"
    inter_binSize = 50000
    intra_binSize = 50000
    win = 10
    reg = 0.2
    cutoff = 0.6
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    result = run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name,
                    output)
    assert result == "File does not exist."


def test_cutoff():
    # error paramter cutoff
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "ref.len"
    inter_binSize = 50000
    intra_binSize = 50000
    win = 10
    reg = 0.2
    cutoff = 2
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    result = run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name,
                      output)
    assert result == "Parameter Error."


def test_win():
    # error parameter win
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "ref.len"
    inter_binSize = 50000
    intra_binSize = 50000
    win = 11.5
    reg = 0.2
    cutoff = 0.6
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    result = run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name,
                      output)
    assert result == "Parameter Error."


def test_interBinSize():
    # error binSize 
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "ref.len"
    inter_binSize = 100000
    intra_binSize = 50000
    win = 10
    reg = 0.2
    cutoff = 0.6
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    result = run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name,
                      output)
    assert result == "Matrix file Error."


def test_matrixConsisRef():
    ## The reference file does not match the matrix file
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "hg38.len"
    inter_binSize = 50000
    intra_binSize = 50000
    win = 10
    reg = 0.2
    cutoff = 0.6
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    result = run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name,
                      output)
    assert result == "Matrix file Error."


def test_matrixFile():
    ## intra_matrix file not exist
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "ref.len"
    inter_binSize = 50000
    intra_binSize = 50000
    win = 10
    reg = 0.2
    cutoff = 0.6
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix_test"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    result = run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name,
                      output)
    assert result == "Matrix file Error."


def test_NormalSample():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    output = current_directory + "/result"
    ref_file = "ref.len"
    inter_binSize = 50000
    intra_binSize = 50000
    win = 10
    reg = 0.2
    cutoff = 0.6
    n_cores = 20
    intra_hicfile = current_directory + "/Matrix_data/Intra_matrix"
    inter_hicfile = current_directory + "/Matrix_data/Inter_matrix"
    name = "K562"
    run_HiSV(ref_file, inter_binSize, intra_binSize, win, reg, cutoff, n_cores, intra_hicfile, inter_hicfile, name, output)
    intra_GT_file = "K562_intra_GT.txt"
    inter_GT_file = "K562_inter_GT.txt"
    intra_result_file = current_directory + "/result/HiSV_intra_SV_result.txt"
    inter_result_file = current_directory + "/result/HiSV_inter_SV_result.txt"
    intra_GT = read_intra_groundtruth(intra_GT_file)
    inter_GT = read_inter_groundtruth(inter_GT_file)
    intra_f1_score = sta_intra_result(intra_GT, intra_result_file)
    inter_f1_score = sta_inter_result(inter_GT, inter_result_file, binSize=50000)
    result_score = (inter_f1_score + intra_f1_score) / 2
    assert result_score > 0.5


if __name__ == '__main__':
    pytest.main(["-v", "HiSV_test.py"])
