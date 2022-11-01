from optparse import OptionParser
import coolToMatrix
import hicToMatrix
import bamToMatrix
import os

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-o', '--output', default='', help='dir of output file')
    parser.add_option('-m', '--hicfile', default='', help='input hic file')
    parser.add_option('-t', '--type', default='', help='type of input hic file')
    parser.add_option('-r', '--length', default='', help='length of reference genome file')
    parser.add_option('-b', '--binSize', default=50000, type ='int', help='binSize of final matrix file')
    (opts, args) = parser.parse_args()
    output = opts.output
    if not os.path.exists(output):
        os.makedirs(output)
    hicfile = opts.hicfile
    ref_file = opts.length
    name = 'sample'
    binSize = opts.binSize
    covert_type = opts.type
    type_list = ['mcool', 'cool', 'hic', 'bam']
    if covert_type not in type_list:
        print("")
    else:
        if covert_type == 'mcool':
            hicfile = hicfile + '::/resolutions/' + str(binSize)
            print(hicfile)
            coolToMatrix.coolToMatrix(hicfile, binSize, output, name)
        elif covert_type == 'cool':
            coolToMatrix.coolToMatrix(hicfile, binSize, output, name)
        elif covert_type == 'hic':
            hicToMatrix.hicToMatrix(hicfile, binSize, ref_file, output, name)
        else:
            bamToMatrix.bamToMatrix(hicfile, ref_file, binSize, output, name)
