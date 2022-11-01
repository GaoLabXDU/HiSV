import os,sys
import numpy as np
from scipy.sparse import coo_matrix
import straw


def get_chromInfo(chromlf):
    chroms = []
    infos = {}
    inf = open(chromlf)
    for line in inf:
        line = line.strip().split('\t')
        infos[line[0]] = int(line[1])
        chroms.append(line[0])
    return chroms,infos


def dumpMatrix(chrom1, chrom2, binSize, hicfile, chromInfo, outputMatrixFile):
    chrom1length = chromInfo[chrom1]
    chrom2length = chromInfo[chrom2]
    binnumber1 = int(np.divide(chrom1length,binSize)) + 1
    binnumber2 = int(np.divide(chrom2length,binSize)) + 1
    chr1 = chrom1.lstrip('chr')
    chr2 = chrom2.lstrip('chr')
    print(chr1, chr2)
    # result = straw.straw('NONE', '/home/li/NAS/Hi-C_Data/Cell_lines/k562/GSE63525_K562_combined.hic', 'X', 'X', 'BP',1000000)
    result = straw.straw('NONE', hicfile, str(chr1), str(chr2), 'BP', binSize)

    row = np.divide(result[0],binSize)
    col = np.divide(result[1],binSize)
    data = result[2]
    # print(max(row),max(col),binnumber1,binnumber2)
    #print max(row),max(col),binnumber1,binnumber2
    res = coo_matrix((data, (row, col)), shape=(binnumber1, binnumber2)).toarray()
    if chrom1 == chrom2:
        res += res.T - np.diag(res.diagonal())
    np.savetxt(outputMatrixFile, res, fmt='%i', delimiter='\t')


def hicToMatrix(hicfile, binSize, ref_file, outdir,name):
    MatrixInfo = {}
    chroms, chromInfo = get_chromInfo(ref_file)
    chr_count = len(chroms)
    outdir_inter = os.path.join(outdir, "Inter_matrix")
    outdir_intra = os.path.join(outdir, "Intra_matrix")
    if not os.path.isdir(outdir_inter):
        os.mkdir(outdir_inter)
    if not os.path.isdir(outdir_intra):
        os.mkdir(outdir_intra)

    for i in range(3, 4):
        for j in range(i, i+1):
            if i == j:
                outdir = outdir_intra
            else:
                outdir = outdir_inter
            chrom1 = chroms[i]
            chrom2 = chroms[j]
            outputname = os.path.join(outdir,name + '_%skb_%s_%s_matrix.txt' % (str(binSize//1000), chrom1, chrom2))
            dumpMatrix(chrom1, chrom2, binSize, hicfile, chromInfo, outputname)
            MatrixInfo[chrom1 + '_' + chrom2] = outputname
    return MatrixInfo


ref_file = "/home/li/PycharmProjects/HiSV/hg19.len"
hicfile = "/home/li/NAS/Hi-C_Data/Cell_lines/k562/GSE63525_K562_combined.hic"
binSize = 50000
outdir = "/media/li/Data/K562/TAD_result"
name = "k562"
hicToMatrix(hicfile, binSize, ref_file, outdir, name)
