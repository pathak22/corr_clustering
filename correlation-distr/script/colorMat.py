import pylab
from numpy import array, zeros

import sys

from evaluate import readClMat, readDataFile

def cluster(inds):
    mat = zeros((len(inds), len(inds)))

    for ii,cli in enumerate(inds):
        mat[ii][ii] = 1
        for jj,clj in enumerate(inds):
            if cli.label == clj.label:
                mat[ii][jj] = 1
                mat[jj][ii] = 1
            if jj > ii:
                break
    return mat

if __name__ == "__main__":
    try:
        mat = readClMat(file(sys.argv[1]))

        if len(mat) == 1:
            raise ValueError("Probably a data file")

    except ValueError:
        dataFh = file(sys.argv[1])
        data = readDataFile(dataFh)
        data2 = readDataFile(dataFh)
        if data2:
            mat = cluster(data2)
        else:
            mat = cluster(data)
    #print mat
    mat = array(mat)
#    print mat
    pylab.axis([0, len(mat), len(mat), 0])
    #    pylab.gray()
    pylab.pcolor(mat)
    pylab.colorbar()
    pylab.show()
