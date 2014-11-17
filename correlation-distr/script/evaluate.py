from __future__ import division

import sys
from path import path

from Classifiers import Classifier, readDataFile

from math import sqrt, log
from AIMA import DefaultDict
from StringIO import StringIO

from ClusterMetrics import ConfusionMatrix
from Probably import entropy_of_multinomial

def readClMat(ff):
    nodes = int(ff.next())

    mat = []
    for x in range(nodes):
        newRow = [0 for y in range(nodes)]
        mat.append(newRow)

    for row in range(nodes):
        mat[row][row] = 1
        for col in range(row + 1, nodes):
            nextIn = float(ff.next())
            mat[row][col] = nextIn
            mat[col][row] = nextIn

    return mat

def objective(points, mat, log=False):
    """The total weight of cross edges plus the total weight of in-cluster
    anti-edges. Larger is worse."""
    obj = 0

    if log:
        weightFn = safeLog
    else:
        weightFn = lambda x: x

    for same,p1,p2 in Classifier.instances(points):
        weight = mat[p1.index][p2.index]

        if same:
            obj += weightFn(1 - weight)
        else:
            obj += weightFn(weight)

    return obj

def safeLog(weight):
    epsilon = 1e-5
    if weight < epsilon:
        weight = epsilon

    return log(weight)

def logObjective(points, mat):
    """The weight of log cross edges."""
    return objective(points, mat, log=True)

def trivialUpperBound(mat, log=False):
    """Objective value of upper bound of all singletons."""
    if log:
        weightFn = safeLog
    else:
        weightFn = lambda x: x

    obj = 0

    for i in range(len(mat)):
        for j in range(i + 1, len(mat)):
            obj += weightFn(mat[i][j])
    return obj

def trivialLowerBound(mat, log=False):
    """Objective value of sum over edges of min of both weights."""
    obj = 0

    if log:
        weightFn = safeLog
    else:
        weightFn = lambda x: x

    for i in range(len(mat)):
        for j in range(i + 1, len(mat)):
            obj += min(weightFn(mat[i][j]), weightFn(1 - mat[i][j]))
    return obj

class Eval:
    def __init__(self, gold=None, prop=None, clMat=None, filename="unknown"):
        self.filename = filename
        self.stats = {}

        if gold != None:
            self.stats["clusters"] = len(set([pt.label for pt in prop]))
            self.stats["objective"] = objective(prop, clMat)
            self.stats["logObjective"] = logObjective(prop, clMat)

            conf = ConfusionMatrix()

            for pt,propPt in zip(gold, prop):
                conf.add(pt.label, propPt.label)

            self.stats["rand"] = conf.rand_index()
            self.stats["jaccard"] = conf.jaccard_index()
            self.stats["mirkin"] = conf.mirkin_metric()
            (prec,rec,f) = conf.prec_rec()
            self.stats["prec"] = prec
            self.stats["rec"] = rec
            self.stats["f"] = f
            self.stats["121"] = conf.one_to_one_optimal(verbose=False)[-1]
            self.stats["m21"] = conf.many_to_one(verbose=False)[-1]
            self.stats["vi"] = conf.variation_of_information()
            self.stats["vibound"] = conf.variation_of_information_upper_bound()
            self.stats["nmi"] = conf.normalized_mutual_information()

    def __repr__(self):
        return "File: %s" % self.filename + \
"""
Clusters: %(clusters).3g
Objective: %(objective).3g
Objective (log): %(logObjective).3g

Some edge-counting metrics:
Rand index (max 1): %(rand)g
Jaccard index (max 1): %(jaccard)g
Mirkin metric (min 0): %(mirkin)g
(Same cluster) Prec: %(prec).3g Rec: %(rec).3g F: %(f).3g

Some node-counting metrics:
One-to-one match (max 1): %(121).3g
Many-to-one match (max 1): %(m21).3g
Variation of information (min 0, max %(vibound).3g): %(vi).3g
Normalized mutual information (0-1): %(nmi).3g
""" % self.stats + \
    self.chatString()

    def chatString(self):
        if "chatScores" not in self.__dict__:
            return ""
        else:
            res = StringIO()
            for metric,mscores in self.chatScores.items():
                print >>res, metric
                print >>res, "\tmin", min(mscores)
                print >>res, "\tmax", max(mscores)
                print >>res, "\tmean", sum(mscores)/len(mscores)
            return res.getvalue()

    def __add__(self, other):
        res = Eval()
        for key,val in self.stats.items():
            res.stats[key] = val + other.stats[key]
        return res

    def normalized(self, norm):
        res = Eval()
        for key,val in self.stats.items():
            res.stats[key] = val / norm
        return res

if __name__ == "__main__":
    try:
        gold = path(sys.argv[1])
        classFileName = path(sys.argv[2])

        goldFile = file(gold)
        classFile = file(classFileName)
    except IndexError:
        print "eval [gold node labels] [classifier output] " +\
              "[induced node labels]"
        sys.exit(1)

    train = readDataFile(goldFile)
    test = readDataFile(goldFile)

    clMat = readClMat(classFile)

    print "True clustering has", \
          len(set([pt.label for pt in test])), "clusters"
    print "Objective value of truth:", objective(test, clMat)

    props = []
    for propFile in sys.argv[3:]:
        prop = readDataFile(file(propFile))

        assert(len(prop) == len(test) == len(clMat))

        ev = Eval(test, prop, clMat, filename=propFile)
        props.append(ev)

    print "Best Rand:"
    print max(props, key=lambda x: x.stats["rand"])

    print "Best Objective:"
    print min(props, key=lambda x: x.stats["objective"])

    print "Average:"
    print reduce(lambda x,y: x+y, props).normalized(len(props))
