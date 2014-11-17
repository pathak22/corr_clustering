from __future__ import division
import optparse
import sys

from math import log, exp
from numpy.random.mtrand import dirichlet, multinomial

def labels(num, clusterSizes):
    cluster = 0
    jj = 0
    for ii in range(num):
        jj += 1
        while jj > clusterSizes[cluster]:
            jj = 1
            cluster = cluster + 1

        yield cluster

def addOptions(parseopts):
    parseopts.add_option("-n", "--number", dest="num", default=1,
                         type="int",
                         help="number of points to generate")
    parseopts.add_option("-t", "--train", dest="trainNum", default=1,
                         type="int",
                         help="number of training points to generate")
    parseopts.add_option("-d", "--dirichlet", dest="clusterAlpha", default=1,
                         type="float",
                         help="prior for cluster size distribution")
    parseopts.add_option("-b", "--balanced", dest="balanced", default=False,
                         action="store_true",
                         help="generate balanced cluster sizes")
    parseopts.add_option("-k", "--k-clusters", dest="numClusters",
                         default=None,
                         type="int",
                         help="number of clusters (default is log(N, 1.5))")

def computePrior(options):
    num = options.num
    clusterAlpha = options.clusterAlpha
    balanced = options.balanced
    numClusters = options.numClusters

    assert(num > 0)

    if numClusters is None:
        numClusters = max(int(log(num, 1.5)), 2)

    print >>sys.stderr, numClusters, "clusters"

    clusterPrior = [clusterAlpha] * numClusters
    if not balanced:
        clusterPrior = dirichlet(clusterPrior)
    else:
        norm = sum(clusterPrior)
        clusterPrior = [xx/norm for xx in clusterPrior]

    return clusterPrior


if __name__ == "__main__":
    parseopts = optparse.OptionParser()
    addOptions(parseopts)

    (options, args) = parseopts.parse_args()
    num = options.num
    trainNum = options.trainNum

    clusterPrior = computePrior(options)

    print >>sys.stderr, "cluster prior:",\
          " ".join(["%.3g" % xx for xx in clusterPrior])

    clusterSizes = multinomial(trainNum, clusterPrior)

    print >>sys.stderr, "training cluster sizes:",\
          " ".join([str(xx) for xx in clusterSizes])

    for label in labels(trainNum, clusterSizes):
        print label

    print

    clusterSizes = multinomial(num, clusterPrior)

    print >>sys.stderr, "cluster sizes:",\
          " ".join([str(xx) for xx in clusterSizes])

    for label in labels(num, clusterSizes):
        print label
