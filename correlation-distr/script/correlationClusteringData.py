from sys import argv
from itertools import izip, count

from analysis.chatStats import *
import sys

if __name__ == "__main__":
    chat = readChat(argv[1])
    predictions = file(argv[2])
    keys = file(argv[3])

    outLabel = file(argv[4], 'w')
    outMat = file(argv[5], 'w')

    sparseMat = {}
    lineMapping = {}

    print >>outLabel, "00"
    print >>outLabel

    ct = count()
    for lineNum,line in enumerate(chat):
        if line.thread != -1:
            lineMapping[lineNum] = ct.next()
            print >>outLabel, line.thread

    nodes = ct.next()

    for predictionLine, keyLine, in izip(predictions, keys):
        prob = float(predictionLine.split()[1])
        (i, j) = map(int, keyLine.split())

        j = lineMapping[j]
        i = lineMapping[i]            

        sparseMat[(i,j)] = prob

    print >>outMat, nodes
    for ii in range(nodes):
        for jj in range(ii + 1, nodes):
            if (jj, ii) in sparseMat:
                print >>outMat, sparseMat[(jj, ii)]
            else:
                print >>outMat, .5
