from __future__ import division
from megam import Megam
from Classifiers import Classifier

import cPickle as pickle

from newsgroup import NewsgroupPost

from AIMA import DefaultDict

from random import random
from math import log, exp, sqrt
from path import path
import sys

class MaxEntNewsgroups(Classifier):
    """This classifier runs logistic regression on the tf-df-bucketed
    coordinate differences (requires the word frequency table to be
    externally loaded)."""
    def __init__(self, train, label, extra):
        Classifier.__init__(self, train, label, extra)
        nfile = path(extra.newsgroups)
        print >>sys.stderr, "dump file:", nfile
        ff = file(nfile)
        self.voc = pickle.load(ff)
        self.revVoc = dict([(v,k) for k,v in self.voc.items()])
        gps = pickle.load(ff)
        self.wfr = pickle.load(ff)
        self.dfs = pickle.load(ff)

        self.classifier = Megam(mode="binomial")

    def cosine(self, p1, p2):
        dot = 0
        norm1 = 0
        norm2 = 0

        for feat, val in p1.feats.items():
            if feat.startswith("LSA"):
                norm1 += val ** 2

                val2 = p2.feats.get(feat, 0)
                dot += val * val2

        for feat, val in p2.feats.items():
            if feat.startswith("LSA"):
                norm2 += val ** 2

        return dot / (sqrt(norm1) * sqrt(norm2))

    def makeFeatDict(self, p1, p2, verbose=False):
        shared = DefaultDict([])
        unshared = DefaultDict([])
        sharedSubj = DefaultDict([])
        unsharedSubj = DefaultDict([])

        #df-buckets
        for feat in set(p1.feats.keys() + p2.feats.keys()):
            if feat.startswith("LSA"):
                continue

            if feat.startswith("SUBJ"):
                fnum = int(feat.lstrip("SUBJ"))
                bucket = int(log(self.dfs[fnum], 1.8))
                val1 = p1.feats.get(feat, 0)
                val2 = p2.feats.get(feat, 0)

                if verbose and val1 and val2:
                    print >>sys.stderr, "SUBJ", self.revVoc[fnum], bucket,

                if (val1 or val2) and not (val1 and val2):
                    unsharedSubj[bucket].append(max(val1, val2)//2)
                else:
                    sharedSubj[bucket].append(min(val1, val2)//2)
                
                continue

            fnum = int(feat)

            if fnum not in self.dfs:
                continue

            if self.wfr[fnum] < 3:
                continue

            bucket = int(log(self.dfs[fnum], 1.8))

            val1 = p1.feats.get(feat, 0)
            val2 = p2.feats.get(feat, 0)

            if verbose and val1 and val2:
                print >>sys.stderr, self.revVoc[fnum], bucket,

            if (val1 or val2) and not (val1 and val2):
                unshared[bucket].append(max(val1, val2)//2)
            else:
                shared[bucket].append(min(val1, val2)//2)

        if 0:
            #tf bucketing
            fdict = DefaultDict(0)
            for bucket,ct in shared.items():
                for freq in ct:
                    if freq > 10:
                        freq = 10
                    fdict["SHARE_%g_%g" % (bucket, freq)] += 1
            for bucket,ct in unshared.items():
                for freq in ct:
                    if freq > 10:
                        freq = 10
                    fdict["UNIQUE_%g_%g" % (bucket, freq)] += 1
        elif 0:
            #only df bucketing
            fdict = DefaultDict(0)
            for bucket,ct in shared.items():
                items = len(ct)
                fdict["SHARE_%g" % bucket] += items
            for bucket,ct in unshared.items():
                items = len(ct)
                fdict["UNIQUE_%g" % bucket] += items
        else:
            #proportion features
            fdict = DefaultDict(0)
            for bucket,ct in shared.items():
                items = len(ct)
                nUnshared = len(unshared[bucket])
                frac = items/(nUnshared + items)

                fdict["PROP_%g" % bucket] = frac
                if nUnshared + items == 0:
                    fdict["NO_%g" % bucket] = 1
            for bucket,ct in sharedSubj.items():
                items = len(ct)
                nUnshared = len(unsharedSubj[bucket])
                frac = items/(nUnshared + items)

                fdict["PROP_SUBJ_%g" % bucket] = frac
                if nUnshared + items == 0:
                    fdict["NO_SUBJ_%g" % bucket] = 1

        #LSA
        cos = self.cosine(p1, p2)
        if verbose:
            print >>sys.stderr, "cosine", cos
        fdict["COS"] = cos

        if verbose and 0:
            print >>sys.stderr
            print >>sys.stderr, "Feats", fdict

        if verbose:
            print >>sys.stderr

        return fdict

    def addTrainingInst(self, instance):
        (label, p1, p2) = instance
        fdict = self.makeFeatDict(p1, p2)
        if label == 0:
            label = .2
        self.classifier.add(label, fdict)

    def estimate(self):
        self.classifier.train()

    def classify(self, instance, verbose=False):
        (label, p1, p2) = instance
        if verbose:
            print >>sys.stderr, label, p1.label, p2.label
        fdict = self.makeFeatDict(p1, p2, verbose=verbose)

        res = self.classifier.classify(fdict)
        res = res
        gold = label

        if verbose:
            print >>sys.stderr, res

        self.scores.add(int(gold > .5), int(res > .5))
        return res

if __name__ == "__main__":
    nfile = path("/u/melsner/research/correlation/20ng/groups.dump")
    print >>sys.stderr, nfile
    ff = file(nfile)
    voc = pickle.load(ff)
    gps = pickle.load(ff)
    wfr = pickle.load(ff)
    dfs = pickle.load(ff)

    class WW:
        def __init__(self):
            self.newsgroups = str(nfile)
    extra = WW()
    x = MaxEntNewsgroups([], [], extra)

    x.train()
    print x
    pickle.dump(x, file("/dev/null", 'w'))
