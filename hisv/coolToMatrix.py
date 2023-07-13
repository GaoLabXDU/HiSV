import os, sys
import numpy as np
import cooler
import multiprocessing


def getBins(coolfile):
    """
     Read the header file of the cool file
    :param coolfile: input cool file
    :return: binsInfor {'chrom'：[start_pos, end_pos]}
    """

    binsInfo = {}
    chroms = coolfile.chroms()["name"][:]
    for chrom in chroms:
        idxarray = np.where(coolfile.bins()["chrom"][:] == chrom)
        chromstart = idxarray[0][0]
        chromend = idxarray[0][-1]
        binsInfo[chrom] = [chromstart, chromend]
    return binsInfo


def get_chromInfo(chromlf):

    """
    read reference genome length file
    :param chromlf: reference genome length file
    :return: chromID and chromosome length
    """
    chroms = []
    inf = open(chromlf)
    for line in inf:
        line = line.strip().split('\t')
        chroms.append(line[0])
    return chroms

def dumpMatrix(binSize, coolfile, binsInfo, chrom1, chrom2, outdir, name):
    """
    Convert cool file to contact matrix for each chromosome
    :param binSize: the resolution of Hi-C contact matrix
    :param coolfile: input cool file
    :param binsInfo: chromosome info in cool file
    :param chrom1: chrom1 id
    :param chrom2: chrom2 id
    :param outdir: path to the output file
    :param name: sample name
    :return: current Hi-C matrix
    """

    chrom1start, chrom1end = binsInfo[chrom1]
    chrom2start, chrom2end = binsInfo[chrom2]
    if 'chr' not in chrom1:
        outputname = os.path.join(outdir, name + '_%skb_%s_%s_matrix.txt' % (str(binSize // 1000),
                                                                             'chr'+chrom1, 'chr'+chrom2))
    else:
        outputname = os.path.join(outdir, name + '_%skb_%s_%s_matrix.txt' % (str(binSize // 1000), chrom1, chrom2))

    # matrixName = '_'.join([name, str(resolution)+'kb',chrom1,chrom2,"InterMap_matrix.txt"])
    # matrixfile = os.path.join(outdir,matrixName)
    # print(chrom1start, chrom1end, chrom2start, chrom2end)
    matrix = coolfile.matrix(balance=False)[chrom1start:(chrom1end + 1), chrom2start:(chrom2end + 1)]
    np.savetxt(outputname, matrix, fmt='%.2f', delimiter='\t')


def coolToMatrix(matrixFile, ref_file, binSize, outdir, name, n_cores):

    """
    :param matrixFile: input cool file
    :param binSize: the resolution of Hi-C contact matrix
    :param outdir: path to the output file
    :param name: sample name
    """

    coolfile = cooler.Cooler(matrixFile)
    binsInfo = getBins(coolfile)
    chr_id = get_chromInfo(ref_file)

    chr_count = len(chr_id)
    outdir_inter = os.path.join(outdir, "Inter_matrix")
    outdir_intra = os.path.join(outdir, "Intra_matrix")
    if not os.path.isdir(outdir_inter):
        os.mkdir(outdir_inter)
    if not os.path.isdir(outdir_intra):
        os.mkdir(outdir_intra)
    p = multiprocessing.Pool(processes=n_cores)

    for i in range(chr_count):
        for j in range(i, chr_count):
            if i == j:
                outdir = outdir_intra
            else:
                outdir = outdir_inter

            if 'chr' not in next(iter(binsInfo)):
                chrom1 = chr_id[i].lstrip('chr')
                chrom2 = chr_id[j].lstrip('chr')
            else:
                chrom1 = chr_id[i]
                chrom2 = chr_id[j]
            print('chroms: ', chrom1, chrom2)
            # dumpMatrix(binSize, coolfile, binsInfo, chrom1, chrom2, outdir, name)
            p.apply_async(dumpMatrix, args=(binSize, coolfile, binsInfo, chrom1, chrom2, outdir, name))
    p.close()
    p.join()

