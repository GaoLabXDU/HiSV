import numpy as np
import argparse, sys, os

def read_reflen(chr_len_file, binSize):
    '''
    read chrom length file
    :param chr_len_file:
    :return: chrom ID, chrom bin count
    '''
    with open(chr_len_file) as f:
        chr_id = []
        chr_bin_num = []
        for line in f:
            line_data = line.strip().split('\t')
            chr_id.append(line_data[0])
            chr_bin_num.append(int(int(line_data[1]) / binSize))

    chr_bin_num = np.array(chr_bin_num)
    return chr_id, chr_bin_num


def remove_gap_region(matrix, region_start, region_end):
    '''

    :param intra-chromosomal contact matrix:
    :param gap region strat position:
    :param gap region end position:
    :return:
    '''
    win = 10
    for i in range(len(region_start)):
        start_pos = region_start[i] - win
        end_pos = region_end[i] + win
        matrix[start_pos:end_pos,] = 0
        matrix[:, start_pos:end_pos] = 0
    return matrix


def read_gap_file(gap_file, cur_chr, binSize):
    '''
    :param gap_file:
    :param chrom id:
    :return:
    '''
    region_start = []
    region_end = []
    with open(gap_file) as f:
        for line in f.readlines():
            linelist = line.strip().split('\t')
            if linelist[1] == cur_chr:
                region_start.append(int(linelist[2])//binSize)
                region_end.append(int(linelist[3])//binSize)
    return region_start, region_end


def getargs():
    ## Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(description='''Remove the effect of gap region in hg19 on contact maps.''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--hic_file', help='Path of hi-c contact matrix.')
    parser.add_argument('--binSize', type=int, default=50000, help='Resolution of Hi-C contact matrix.')
    parser.add_argument('--gap_file', help='Gap region file of hg19.')
    parser.add_argument('--ref', help='Reference Genome length file')
    parser.add_argument('--output', help='''Output file path.''')
    parser.add_argument('--name', help='Sample name')


    ## Parse the command-line arguments
    commands = sys.argv[1:]
    if not commands:
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands


def run():
    args, commands = getargs()
    output = args.output
    if not os.path.exists(output):
        os.makedirs(output)
    hic_file = args.hic_file
    binSize = args.binSize
    gap_file = args.gap_file
    chr_len_file = args.ref
    name = args.name
    print("Read chrom length file...")
    chr_id, chr_bin_num = read_reflen(chr_len_file, binSize)
    print(chr_id, chr_bin_num)
    # chr_count = len(chr_id)

    for i in range(len(chr_id)):

        matrix_file = os.path.join(hic_file, name + '_%skb_%s_%s_matrix.txt' % (str(binSize // 1000),
                                                                                chr_id[i], chr_id[i]))
        print("Load matrix file: ", matrix_file)
        contact_matrix = np.loadtxt(matrix_file)
        # raw_num = col_num = np.shape(contact_matrix)[0]
        region_start, region_end = read_gap_file(gap_file, chr_id[i], binSize)
        # print(region_start, region_end)
        remove_matrix = remove_gap_region(contact_matrix, region_start, region_end)
        result_file = os.path.join(output, name + '_%skb_%s_%s_matrix.txt' % (str(binSize // 1000),
                                                                                        chr_id[i], chr_id[i]))

        np.savetxt(result_file, remove_matrix, fmt='%i')


if __name__ == '__main__':
    run()

