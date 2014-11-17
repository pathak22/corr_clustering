from __future__ import division

from Hogwash import Session
from Hogwash.Results import ResultsFile
from waterworks.Processes import bettersystem
import pylab

import sys
from path import path
import os
from StringIO import StringIO
import re
from AIMA import DefaultDict
    
def colorType(model):
    modelStr = " ".join(model)

    if "lp" in modelStr:
        line = "--"
    elif "sdp" in modelStr:
        line = ":"
    else:
        line = "-"

    if "pivot" in modelStr:
        color = "r"
    elif "first" in modelStr:
        color = "c"
    elif "best" in modelStr:
        color = "g"
    elif "vote" in modelStr:
        color = "b"
    else:
        color = "k"

    if "boem" in modelStr:
        sym = "."
    else:
        sym = ""

    return color+sym+line

if __name__ == "__main__":
    workdir = path(sys.argv[1])
    print "Working directory", workdir
    ssnName = "hog" + workdir.dirname().basename()
    print "Session", ssnName
    session = Session(ssnName, read_only=True)
    evalSsnName = "hogeval" + workdir.dirname().basename()
    print "Eval session", evalSsnName
    evalSsn = Session(evalSsnName, read_only=True)

    plotStat = "time"

    modelSeqs = DefaultDict([])

    for ct,job in enumerate(evalSsn):
        if ct % 20 == 0:
            print ct, "..."
        
        if job.status == "finished":
            res = job.results

            try:
                time = res.stats["time"]
            except KeyError:
                time = res.time
            if time == None:
                continue
            
            dset = res.filename.basename()
            jobSize = re.match("^data(\d+)-(\d+)$", dset).group(1)
            jobSize = int(jobSize)

            try:
                jobVal = res.stats[plotStat]
            except KeyError:
                jobVal = res.__dict__[plotStat]

            model = tuple(job.args[0].args[0])

            modelSeqs[model].append((jobSize, jobVal))

    legendItems = []

    print "aggregating"

    for model,seq in modelSeqs.items():
        aggregate = DefaultDict([])
        for size,val in seq:
            aggregate[size].append(val)

        xs = []
        ys = []
        
        for size,vals in sorted(aggregate.items(), key=lambda x: x[0]):
            if None not in vals:
                newVal = sum(vals)/len(vals)
                xs.append(size)
                ys.append(newVal)

        legendItems.append(model)
        #print xs
        #print ys
        pylab.plot(xs, ys, colorType(model))
    pylab.legend(legendItems)
    pylab.show()
