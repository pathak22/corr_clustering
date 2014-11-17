from path import path
import sys
import re

import cPickle as pickle

from itertools import izip, chain
from iterextras import batch

from newsgroup import NewsgroupPost

def symbolize(pt, names):
    name = names.index(pt.group)

    res = ("%d " % name) + \
          " ".join(["%d %d" % (k,v) for k,v in pt.counts.items()]) + \
          " " + \
          " ".join(["SUBJ%d %d" % (k,v) for k,v in pt.subjCounts.items()])
    return res

if __name__ == "__main__":
    ngDump = path("WORD FREQS DUMP")
    fh = file(ngDump)

    vocab = pickle.load(fh)
    groups = pickle.load(fh)
    wfreqs = pickle.load(fh)

    groupNames = groups.keys()
#     for ct,x in enumerate(groupNames):
#         print ct, x
#     sys.exit(0)

    holdout = 4

    outdir = path("PLACE TO WRITE FILES")

    for ct,tts in enumerate(batch(range(len(groups)), holdout)):
        print [groupNames[ti] for ti in tts]

        outfile = outdir/("ng%d" % ct)
        print outfile
        outfile = file(outfile, 'w')

        for pt in chain(*[groups[groupNames[ti]] for ti in tts]):
            print >>outfile, symbolize(pt, groupNames)

        print >>outfile

        for pt in chain(*[groups[groupNames[ti]]
                         for ti in range(len(groupNames))
                         if ti not in tts]):
            print >>outfile, symbolize(pt, groupNames)
