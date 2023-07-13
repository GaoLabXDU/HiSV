import os,sys
import numpy as np
import hicstraw
import multiprocessing


def get_chromInfo(chromlf):

    """
    read reference genome length file
    :param chromlf: reference genome length file
    :return: chromID and chromosome length
    """
    chroms = []
    infos = {}
    inf = open(chromlf)
    for line in inf:
        line = line.strip().split('\t')
        infos[line[0]] = int(line[1])
        chroms.append(line[0])
    return chroms, infos


def dumpMatrix(chrom1, chrom2, binSize, hicfile, chromInfo, outputMatrixFile):
    """
    convert hic file to matrix file for each chromosome
    :param chrom1: chrom1 ID
    :param chrom2: chrom2 ID
    :param binSize: the resolution of Hi-C contact matrix
    :param hicfile:  input hic file
    :param chromInfo: chromID and length
    :param outputMatrixFile: path to the output file
    """
    chrom1length = chromInfo[chrom1]
    chrom2length = chromInfo[chrom2]
    binnumber1 = int(np.divide(chrom1length, binSize)) + 1
    binnumber2 = int(np.divide(chrom2length, binSize)) + 1
    # Make chromID consistent with the definition of chromID in the hi-c file
    chr1 = chrom1.lstrip('chr')
    chr2 = chrom2.lstrip('chr')
    print("chromosome: ", chr1, chr2)
    # initialization contact matrix according the length of chromosomes
    init_matrix = np.full((binnumber1, binnumber2), 0)
    # result = hicstraw.straw('NONE', hicfile, str(chrom1), str(chrom2), 'BP', binSize)
    result = hicstraw.straw('observed', 'NONE', hicfile, str(chr1), str(chr2), 'BP', binSize)

    for i in range(len(result)):
        # print(result[i].binX, result[i].binY, result[i].counts)
        row = int(result[i].binX / binSize)
        col = int(result[i].binY / binSize)
        init_matrix[row][col] = result[i].counts

    if chrom1 == chrom2:
        init_matrix += init_matrix.T - np.diag(init_matrix.diagonal())
    np.savetxt(outputMatrixFile, init_matrix, fmt='%i', delimiter='\t')


def hicToMatrix(hicfile, binSize, ref_file, outdir, name, n_cores):
    """
    :param hicfile: input hic file
    :param binSize: the resolution of Hi-C contact matrix
    :param ref_file: reference genome length file
    :param outdir: path to the output file
    :param name: sample name
    :param n_cores: cores number
    """
    MatrixInfo = {}
    chroms, chromInfo = get_chromInfo(ref_file)
    chr_count = len(chroms)
    # Intra-chromosomal and inter-chromosomal contact matrices are stored separately
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

            chrom1 = chroms[i]
            chrom2 = chroms[j]
            outputname = os.path.join(outdir, name + '_%skb_%s_%s_matrix.txt' % (str(binSize//1000), chrom1, chrom2))
            # dumpMatrix(chrom1, chrom2, binSize, hicfile, chromInfo, outputname)
            p.apply_async(dumpMatrix, args=(chrom1, chrom2, binSize, hicfile, chromInfo, outputname))
    p.close()
    p.join()



