from __future__ import division
from megam import Megam
from classify import readDataFile

import cPickle as pickle

from newsgroup import NewsgroupPost

from itertools import count

from random import random
from math import log
from path import path
import sys

if __name__ == "__main__":
    tfile = file(sys.argv[1])

    tpoints = readDataFile(tfile)

    ct = 0

    for point in tpoints:
        for feat,val in point.feats.items():
            print ct, int(feat), val
        ct += 1

    tspoints = readDataFile(tfile)
    for point in tspoints:
        for feat,val in point.feats.items():
            print ct, int(feat), val
        ct += 1

