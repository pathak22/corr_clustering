from __future__ import division
from megam import Megam
from classify import readDataFile

import cPickle as pickle

from newsgroup import NewsgroupPost

from AIMA import DefaultDict

from itertools import count, izip

from random import random
from math import log
from path import path
import sys

if __name__ == "__main__":
    tfile = file(sys.argv[1])
    ufile = file(sys.argv[2])

    tpoints = readDataFile(tfile)
    ct = 0

    for point in tpoints:
        print point.label,
        for ind,coord in enumerate(ufile.next().split()):
            print "LSA%d" % ind, coord,
        for feat,val in point.feats.items():
            print feat, val,
        print
        
        ct += 1

    print

    tspoints = readDataFile(tfile)
    for point in tspoints:
        print point.label,
        for ind,coord in enumerate(ufile.next().split()):
            print "LSA%d" % ind, coord,
        for feat,val in point.feats.items():
            print feat, val,
        print
        
        ct += 1
