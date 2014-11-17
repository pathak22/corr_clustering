from __future__ import division
import optparse
import sys
import itertools

from numpy.random.mtrand import normal, multinomial, gamma
from numpy import zeros, array

from genFeatureless import labels, addOptions, computePrior

def invGamma(alpha, beta):
    return 1.0 / gamma(alpha, 1/beta)

def features(label, means, variances, noisemeans, noisevars):
    feats = []
    fctr = itertools.count()
    for mean,var in zip(means[label], variances[label]):
        feats.append(next(fctr))
        feats.append(normal(mean, var))
    for mean,var in zip(noisemeans, noisevars):
        feats.append(next(fctr))
        feats.append(normal(mean, var))
    return feats

if __name__ == "__main__":
    parseopts = optparse.OptionParser()
    addOptions(parseopts)

    parseopts.add_option("-f", "--features", dest="feats", default=1,
                         type="int",
                         help="number of features per point")
    parseopts.add_option("-z", "--noise", dest="noisefeats", default=0,
                         type="int",
                         help="number of noise features per point")
    parseopts.add_option("-s", "--sigma", dest="sigma", default=1,
                         type="float",
                         help="fixed feature variance")
    parseopts.add_option("-v", "--var", dest="cvar", default=30,
                         type="float",
                         help="variance of distribution of cluster means")
    parseopts.add_option("-a", "--alpha", dest="vvar", default=0,
                         type="float",
                         help="alpha parameter to inv-gamma prior on"+
                         " feature variance")

    (options, args) = parseopts.parse_args()
    num = options.num
    trainNum = options.trainNum
    cvar = options.cvar
    vvar = options.vvar

    feats = options.feats
    noiseFeats = options.noisefeats
    sigma = options.sigma

    print(sys.stderr, "%d real features, %d noise features" %\
          (feats, noiseFeats))

    clusterPrior = computePrior(options)

    #cluster means are normal(0, sigma)
    means = zeros((len(clusterPrior),feats))
    for cluster in range(len(clusterPrior)):
        for feat in range(feats):
            means[cluster][feat] = normal(0, cvar)

    variances = zeros((len(clusterPrior), feats))
    if vvar == 0:
        variances += sigma
    else:
        for cluster in range(len(clusterPrior)):
            for feat in range(feats):
                variances[cluster][feat] += invGamma(vvar, 1.0)

    #currently all noise features have mean 0
    noiseMeans = zeros((noiseFeats,))

    #noise variances are all sigma
    noiseVariances = zeros((noiseFeats,))
    if vvar == 0:
        noiseVariances += sigma
    else:
        for feat in range(noiseFeats):
            noiseVariances[feat] += invGamma(vvar, 1.0)

    print (sys.stderr, "cluster prior:",\
          " ".join(["%.3g" % xx for xx in clusterPrior]))

    print (sys.stderr, "means:")
    for row in means:
        print (sys.stderr, " ".join(["%.3g" % xx for xx in row]))

    if vvar == 0:
        print (sys.stderr, "fixed variance %.3g" % sigma)
    else:
        print (sys.stderr, "variances:")
        for row in variances:
            print (sys.stderr, " ".join(["%.3g" % xx for xx in row]))
        print (sys.stderr, "noise variances:")
        print (sys.stderr, " ".join(["%.3g" % xx for xx in noiseVariances]))

    clusterSizes = multinomial(trainNum, clusterPrior)

    print (sys.stderr, "training cluster sizes:",\
          " ".join([str(xx) for xx in clusterSizes]))

    for label in labels(trainNum, clusterSizes):
        print (label, " ".join([str(xx) for xx in
                               features(label, means, variances,
                                        noiseMeans, noiseVariances)]))

    clusterSizes = multinomial(num, clusterPrior)
    print(clusterSizes)
    print(sys.stderr, "cluster sizes:",\
          " ".join([str(xx) for xx in clusterSizes]))

    for label in labels(num, clusterSizes):
        print(label, " ".join([str(xx) for xx in
                               features(label, means, variances,
                                        noiseMeans, noiseVariances)]))
