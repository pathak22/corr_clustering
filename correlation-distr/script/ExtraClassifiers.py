from __future__ import division

import optparse
import sys

from iterextras import batch
import itertools
from numpy import array, mean, std
from numpy.random.mtrand import beta, random_sample as rand
from scipy.stats import norm
from Classifiers import Classifier

from scores import Scores

def bernoulli(p):
    return rand() < p

class RandomErrClassifier(Classifier):
    """This classifier makes random errors (without regard to the actual
    datapoints involved) on an epsilon fraction of edges.
    Correct in-class link strengths are iid ~ Beta(a-in, b-in),
    and correct cross-class links are iid ~ Beta(a-out, b-out)."""
    def __init__(self, train, label, extra):
        Classifier.__init__(self, train, label, extra)
        self.aIn = extra.aIn
        self.aOut = extra.aOut
        self.bIn = extra.bIn
        self.bOut = extra.bOut
        self.epsilon = extra.epsilon
        self.simple = extra.simple
    def addTrainingInst(self, instance):
        pass

    def estimate(self):
        pass

    def classify(self, instance):
        (label, p1, p2) = instance

        useLabel = label
        if bernoulli(self.epsilon):
            useLabel = not label

        if self.simple:
            res = float(useLabel)
        else:
            if useLabel:
                res = beta(self.aIn, self.bIn)
            else:
                res = beta(self.aOut, self.bOut)
                
        self.scores.add(label, int(res > .5))

        return res

class NaiveBayesClassifier(Classifier):
    """This classifier assumes that feature attributes are continuous,
    and models the difference of each feature as a sample from a Gaussian."""
    def __init__(self, train, label, extra):
        Classifier.__init__(self, train, label, extra)
        self.data = { 0 : {}, 1 : {} }

    def addTrainingInst(self, instance):
        (label, p1, p2) = instance
        for feat in set(p1.feats.keys() + p2.feats.keys()):
            if feat not in self.data[label]:
                self.data[label][feat] = []

            val1 = p1.feats.get(feat, 0)
            val2 = p2.feats.get(feat, 0)

            diff = abs(val1 - val2)
            self.data[label][feat].append(diff)

    def estimate(self):
        self.means = { 0 : {}, 1 : {} }
        self.stdevs = { 0 : {}, 1 : {} }
        for label,tab in self.data.items():
            for feat,diffs in tab.items():
                self.means[label][feat] = mean(diffs)
                self.stdevs[label][feat] = std(diffs)

        self.data = None

    def classify(self, instance):
        (label, p1, p2) = instance
        ps = [1, 1]
        for feat in set(p1.feats.keys() + p2.feats.keys()):
            val1 = p1.feats.get(feat, 0)
            val2 = p2.feats.get(feat, 0)

            diff = abs(val1 - val2)

            for ii in range(len(ps)):
                ps[ii] *= norm.pdf(diff, self.means[ii][feat],
                                        self.stdevs[ii][feat])

        res = ps[1] / sum(ps)

        self.scores.add(label, int(res > .5))

        return res
