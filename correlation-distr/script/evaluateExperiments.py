from __future__ import division

from Hogwash import Session
from Hogwash.Results import ResultsFile
from waterworks.Processes import bettersystem

import sys
from path import path
import os
from StringIO import StringIO
import re
from AIMA import DefaultDict
from itertools import izip

from evaluate import Eval, readClMat
from classify import readDataFile

def jobTime(job):
    log = job.logfile

    for line in log:
        if "Concert exception caught" in line:
            return None
    
    statusdict = job.statusdict
    assert(job.status == "finished")
    return statusdict["finished"] - statusdict["started"]

def evalTask(job):
    assert(job.status == "finished")

    ddir = job.args[1]
    goldFile = file(ddir/"data")
    trainPts = readDataFile(goldFile)
    truth = readDataFile(goldFile)

    matrixFile = file(ddir/"matrix")
    matrix = readClMat(matrixFile)

    propFile = file(str(job.results))
    prop = readDataFile(propFile)

    res = Eval(gold=truth, prop=prop, clMat=matrix, filename=job.args[1])
    res.stats["time"] = jobTime(job)

    return res

if __name__ == "__main__":
    workdir = path(sys.argv[1])
    print "Working directory", workdir
    wbase = workdir.dirname().basename()
    ssnName = "hog" + wbase
    print "Session", ssnName
    session = Session(ssnName, read_only=True)

    args = list(session)
    evalSsn = Session("evaluateExperiments", "evalTask", args,
                      name="hogeval%s" % wbase)
    depDict = {}
    for evJob,expJob in izip(evalSsn,session):
        depDict[evJob] = expJob
    evalSsn.set_dependencies(depDict)
    
