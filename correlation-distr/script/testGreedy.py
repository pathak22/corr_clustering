from __future__ import division

import sys
from path import path
from StringIO import StringIO

from waterworks.Processes import bettersystem

from classify import readDataFile
from evaluate import Eval, readClMat, objective

if __name__ == "__main__":
    tester = path("bin32/greedy")

    try:
        gold = path(sys.argv[1])
        classFileName = path(sys.argv[2])

        goldFile = file(gold)
        classFile = file(classFileName)
    except IndexError:
        print "testGreedy [gold node labels] [classifier output] "
        sys.exit(1)
        
    train = readDataFile(goldFile)
    test = readDataFile(goldFile)

    clMat = readClMat(classFile)

    print "True clustering has", \
          len(set([pt.label for pt in test])), "clusters"
    print "Objective value of truth: %.3g" % objective(test, clMat)
    print

    for algorithm in ["first", "best", "vote", "pivot"]:
        evals = []

        print "================", algorithm, "================="
        print

        for run in range(10):
            cmd = "%s -a %s %s" % (tester, algorithm, classFileName)
            output = StringIO()
            status = bettersystem(cmd, stdout=output, stderr=StringIO())
            assert(status == 0)

            prop = readDataFile(output.getvalue().split("\n"))

            score = Eval(test, prop, clMat, filename=("run%d" % run))
            evals.append(score)
            print "Run", run
            print score

        print algorithm
        print "Best Objective:"
        print min(evals, key=lambda x: x.stats["objective"])

        print "Average:"
        print reduce(lambda x,y: x+y, evals).normalized(len(evals))
        print
