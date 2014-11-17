from __future__ import division

import optparse
import sys

from iterextras import batch
import itertools

from scores import Scores

def bernoulli(p):
    return rand() < p

def readDataFile(dataFile):
    points = []
    ct = itertools.count()
    for line in dataFile:
        if line.strip():
            points.append(Point(line, index=ct.next()))
        else:
            break
    return points

class Point:
    def __init__(self, line, index=None):
        fields = line.split()
        self.label = fields[0]
        self.feats = dict([(x, float(y)) for (x,y) in batch(fields[1:])])
        self.index = index

    def __repr__(self):
        return "%s %s" % (self.label,
                          " ".join(["%s %s" % (kk,vv)
                                    for (kk,vv) in self.feats.items()]))

class Classifier:
    def __init__(self, train, label, extra):
        self.trainSet = train
        self.labelSet = label
        self.scores = Scores()

    def train(self):
        for instance in self.instances(self.trainSet):
            self.addTrainingInst(instance)
        self.estimate()

    def label(self, out=sys.stdout):
        print >>out, len(self.labelSet)
        for instance in self.instances(self.labelSet):
            print >>out, self.classify(instance)

    def partialLabel(self, rows, out=sys.stdout):
        def rowInstances(points):
            for ii in rows:
                print >>sys.stderr, "row", ii
                for jj in range(ii + 1, len(points)):
                    yield (points[ii].label == points[jj].label,
                           points[ii], points[jj])
        for instance in rowInstances(self.labelSet):
            print >>out, self.classify(instance)            

    def check(self, rounds=1000):
        for rnd in range(rounds):
            ii = randrange(len(self.labelSet))
            jj = randrange(len(self.labelSet))

            instance = (self.labelSet[ii].label == self.labelSet[jj].label,
                        self.labelSet[ii], self.labelSet[jj])

            self.classify(instance, verbose=True)

    @classmethod
    def instances(self, points, verbose=False):
        for ii in range(len(points)):
            if verbose:
                print >>sys.stderr, "row", ii
            for jj in range(ii + 1, len(points)):
                yield (points[ii].label == points[jj].label,
                       points[ii], points[jj])

    def addTrainingInst(self, instance):
        raise NotImplementedError("Can't do this. Please subclass!")

    def estimate(self, instance):
        raise NotImplementedError("Can't do this. Please subclass!")

    def classify(self, instance):
        raise NotImplementedError("Can't do this. Please subclass!")
