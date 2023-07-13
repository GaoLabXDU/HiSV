import os, pyBigWig
import numpy as np
import pandas as pd
from collections import Counter
from pygam import PoissonGAM, s


def get_contacts(matrixfile, chrID, binSize):
    '''
    Extract the coverage of each bin from the Hi-C matrix
    :param matrixfile: Hi-C matrix
    :param chrID: chrom ID
    :param binSize: resolution of the hic matrix
    :return: {'chrom', 'start', 'end', 'Coverage'}
    '''
    matrix = np.loadtxt(matrixfile)
    matrix = matrix[:-1, :-1]
    row, col = np.diag_indices_from(matrix)
    matrix[row, col] = 0
    data = pd.DataFrame(columns=['chrom', 'start', 'end', 'Coverage'])
    for k in range(len(matrix[0])):
        coverage = np.sum(matrix[k])
        start = k * binSize
        end = (k + 1) * binSize
        data = data.append({'chrom': chrID, 'start': start, 'end': end, 'Coverage': coverage}, ignore_index=True)
    return data


def signal_from_bigwig(table, bw_fil, name='GC'):
    '''
    Extract the features (GC and mappability) of each bin from the bigwig file
    :param table: chrom info {'chrom', 'start', 'end', 'Coverage'}
    :param bw_fil: bigwig file
    :param name: GC or mappability
    :return: {'chrom', 'start', 'end', 'Coverage', 'GC', 'mappability'}
    '''
    bw = pyBigWig.open(bw_fil)
    arr = []
    for i in table.index:
        row = table.loc[i]
        # (row['chrom'], row['start'], row['end'])
        v = bw.stats(row['chrom'], row['start'], row['end'], type='mean')[0]
        if v is None:
            v = 0
        arr.append(v)

    table[name] = np.r_[arr]

    return table


def count_REsites(table, npz_fil, res):
    '''
    Extract the REsites (GC and mappability) of each bin from the bigwig file
    :param table: chrom info {'chrom', 'start', 'end', 'Coverage'}
    :param npz_fil: RE file
    :param res: resolution of the hic matrix
    :return: {'chrom', 'start', 'end', 'Coverage', 'GC', 'mappability', 'RE'}
    '''
    RE = np.load(npz_fil)

    RE_by_bin = {}
    for chrom in RE:
        tmp = RE[chrom] // res
        RE_by_bin[chrom] = Counter(tmp)

    arr = []
    for i in table.index:
        row = table.loc[i]
        b_i = row['start'] // res
        chrom = row['chrom']
        arr.append(RE_by_bin[chrom][b_i])

    table['RE'] = np.r_[arr]

    return table


def filterZeros(table, cols=['GC', 'Mappability', 'RE']):
    '''
    Filter out bins with feature 0
    :param table: {'chrom', 'start', 'end', 'Coverage', 'GC', 'mappability', 'RE'}
    :param cols:
    :return: feature == 0 index and filtered tabel
    '''
    mask = (table['GC'] != 0) & (table['Mappability'] != 0) & (table['RE'] != 0) & (table['Coverage'] != 0)
    filtered = table[mask]

    return mask, filtered


def normalized_coverage(matrixfile, gc_file, map_file, cutsize_file, binSize, chrID):
    # Normalized the coverage by possion GAM model
    table = get_contacts(matrixfile, chrID, binSize)
    table = signal_from_bigwig(table, gc_file, name='GC')
    # Load mappability scores
    table = signal_from_bigwig(table, map_file, name='Mappability')
    # Count the number of cut sizes for each bin
    table = count_REsites(table, cutsize_file, binSize)
    # Filter out invalid bins
    mask, filtered = filterZeros(table)
    # normalized coverage by generated addition model
    X = filtered[['GC', 'Mappability', 'RE']].values
    y = filtered['Coverage'].values
    gam = PoissonGAM(s(0) + s(1) + s(2), fit_intercept=True).fit(X, y)
    # gam.gridsearch(X, y)
    # logger.info('Output residuals ...')
    residuals = gam.deviance_residuals(X, y)
    residuals = np.array(residuals - residuals.min())
    # residuals = np.nan_to_num(residuals)
    idx = np.where(mask)[0]
    coverage = np.zeros(table.shape[0])
    coverage[idx] = residuals

    return coverage

'''
gc_file = "/mnt/d/HiSV_test_data/normData/hg19_1kb_GC.bw"
mapscore_file = "/mnt/d/HiSV_test_data/normData/hg19_mappability_100mer.1kb.bw"
cutsites_file = "/mnt/d/HiSV_test_data/normData/hg19.HindIII.npz"
# matrixfile = "/mnt/d/HiSV_test_data/test_data/hic_data/chr9.txt"
matrixfile = "/mnt/d/HiSV_test_data/result/Intra_matrix/K562_50kb_chr9_chr9_matrix.txt"
binSize = 50000
chrID = "chr9"
# table = get_contacts(matrixfile, chrID, binSize)
normed_coverage = normalized_coverage(matrixfile, gc_file, mapscore_file, cutsites_file, binSize, chrID)
'''
