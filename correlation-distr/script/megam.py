from __future__ import division
from path import path
from StringIO import StringIO
from math import exp
import sys

from waterworks.Processes import bettersystem
from waterworks.Files import keepable_tempfile
from tempfile import NamedTemporaryFile

megam = path("PATH-TO-/megam_i686.opt")

class Megam:
    def __init__(self, mode="binary"):
        self.trainData = NamedTemporaryFile()
        self.mode = mode
    def add(self, label, fdict):
        fstr = " ".join(["%s %s" % (k,v) for k,v in fdict.items()])
        self.trainData.write("%g %s\n" % (label, fstr))
    def train(self):
        self.trainData.flush()
        cmd = "%s -fvals %s %s" % (megam, self.mode, self.trainData.name)
        err = StringIO()
        out = StringIO()
        bettersystem(cmd, stdout=out, stderr=err)
        self.trainData.close()
        self.trainData = None

        #print err.getvalue()
        #print out.getvalue()

        self.readModel(out.getvalue())
    def readModel(self, output):
        self.model = {}
        for line in output.split("\n"):
            line = line.strip()
            if not line:
                continue
            (coeff,val) = line.split()
            val = float(val)
            self.model[coeff] = val
    def classify(self, fdict, model=None):
        if not model:
            model = self.model
        dot = model["**BIAS**"]

        for feat,val in fdict.items():
             #features unseen in training are skipped
            if feat in model:
                #print >>sys.stderr, feat, model[feat], val, model[feat]*val
                dot += model[feat] * val

        #print >>sys.stderr, "::", dot

        while True:
            try:
                return 1.0 / (1 + exp(-dot))
            except OverflowError:
                dot /= 2

if __name__ == "__main__":
    cls = Megam()

    for case in [ (1, {"a":1}),
                  (1, {"b":1}),
                  (0, {"c":1, "b":1}),
                  (0, {"c":1}),
                  ]:
        label,fdict = case
        cls.add(label, fdict)
    cls.train()

    print cls.model
    print cls.classify({"a":1, "b":1, "c":1})
    #right answer: .462
