import os, sys
import numpy as np
import cooler


def getBins(coolfile):
    binsInfo = {}
    chr_id = []
    chroms = coolfile.chroms()["name"][:]
    for chrom in chroms:
        idxarray = np.where(coolfile.bins()["chrom"][:] == chrom)
        chromstart = idxarray[0][0]
        chromend = idxarray[0][-1]
        binsInfo[chrom] = [chromstart, chromend]
        chr_id.append(chrom)
    return binsInfo, chr_id


def dumpMatrix(binSize, coolfile, binsInfo, chrom1, chrom2, outdir, name):
    chrom1start, chrom1end = binsInfo[chrom1]
    chrom2start, chrom2end = binsInfo[chrom2]

    outputname = os.path.join(outdir, name + '_%skb_%s_%s_matrix.txt' % (str(binSize // 1000), chrom1, chrom2))

    # matrixName = '_'.join([name, str(resolution)+'kb',chrom1,chrom2,"InterMap_matrix.txt"])
    # matrixfile = os.path.join(outdir,matrixName)
    print(chrom1start, chrom1end, chrom2start, chrom2end)
    matrix = coolfile.matrix(balance=False)[chrom1start:(chrom1end + 1), chrom2start:(chrom2end + 1)]
    np.savetxt(outputname, matrix, fmt='%.2f', delimiter='\t')
    return outputname


def coolToMatrix(matrixFile, binSize, outdir, name):
    MatrixInfo = {}
    coolfile = cooler.Cooler(matrixFile)
    binsInfo, chr_id = getBins(coolfile)
    print(binsInfo)
    print(chr_id)
    chr_count = len(chr_id)
    outdir_inter = os.path.join(outdir, "Inter_matrix")
    outdir_intra = os.path.join(outdir, "Intra_matrix")
    if not os.path.isdir(outdir_inter):
        os.mkdir(outdir_inter)
    if not os.path.isdir(outdir_intra):
        os.mkdir(outdir_intra)

    for i in range(chr_count):
        for j in range(i, chr_count):
            if i == j:
                outdir = outdir_intra
            else:
                outdir = outdir_inter
            chrom1 = chr_id[i]
            chrom2 = chr_id[j]
            print(chrom1, chrom2)
            matrixfile = dumpMatrix(binSize, coolfile, binsInfo, chrom1, chrom2, outdir, name)
            MatrixInfo[chrom1 + '_' + chrom2] = matrixfile
    return MatrixInfo

'''
matrixFile = "/home/li/NAS/Hi-C_Data/Cell_lines/k562/K562.mcool::/resolutions/50000"
binSize = 50000
outdir = "/home/li/PycharmProjects/trans_type_to_matrix/test"
name = "simu"
coolToMatrix(matrixFile, binSize, outdir, name)
'''

