from __future__ import division
from Hogwash import Session

import sys
from path import path
import os
import re

from AIMA import DefaultDict

from DistributedExperiment import Experiment, Evaluation

def getCorpus(obj):
    return int(re.search("\d+", obj.exp.mfile.basename()).group(0))

def solverObj(solvers):
    if "log" in solvers:
        return "logObjective"
    return "objective"

if __name__ == "__main__":
    sessionName = sys.argv[1]

    ssn = Session(sessionName, read_only=True)

    evals = DefaultDict(DefaultDict([]))

    for job in ssn:
        #print job
        if job.status != "finished":
            continue
        obj = job.args[0]
        if isinstance(obj, Evaluation):
            corpus = getCorpus(obj)
            solvers = tuple(obj.exp.solvers)
            evals[corpus][solvers].append(job.results)

    avgBest = DefaultDict([])
    avgAvg = DefaultDict([])

    for corpus,values in evals.items():
        print "CORPUS", corpus, "------------------"
        print
        for solvers,props in values.items():
            print solvers
            print "Best Objective:"
            best = min(props, key=lambda x: x.stats[solverObj(solvers)])
            print best
            avgBest[solvers].append(best)
            print "Average:"
            avg = reduce(lambda x,y: x+y, props).normalized(len(props))
            print avg
            avgAvg[solvers].append(avg)
            #print len(props)
            print
            print

    print "ALL", "-----------------------"
    print
    table = []
    for solvers, props in avgBest.items():
        print solvers
        print "Best Objective:"
        best = reduce(lambda x,y: x+y, props).normalized(len(props))
        print best
        print "Avg"
        props2 = avgAvg[solvers]
        avg = reduce(lambda x,y: x+y, props2).normalized(len(props2))
        print avg
        print
        print

        row = [ " ".join(solvers), best.stats[solverObj(solvers)],
                best.stats["rand"], best.stats["f"], best.stats["121"] ]
        table.append(row)

    table = sorted(table, key=lambda x: x[1])

    newtable = [["", "obj", "rand", "f", "121"] ]
    for row in table:
        print row
        print len(row)
        (name, obj, rand, f, s121) = row
        row2 = [name, "%.3g" % obj, "%.4g" % rand, "%.2g" % f, "%.2g" % s121]
        print row2
        newtable.append(row2)

    from TeXTable import texify
    print texify(newtable)

    if 0:
        #correlation analysis
        from numpy import corrcoef
        from scipy.stats import kendalltau, spearmanr

        for stat in ["121", "f", "rand"]:
            ccs = []
            ccs10 = []

            for corpus,values in evals.items():
                allInCorpus = []
                for solvers,props in values.items():
                    stats = [prop.stats[stat] for prop in props]
                    objs = [prop.stats[solverObj(solvers)] for prop in props]
                    paired = zip(stats, objs)
                    allInCorpus += paired

                allInCorpus = sorted(allInCorpus, key=lambda x:x[1]) #by obj
                tenthPerc = allInCorpus[:(len(allInCorpus)//10)]
                #print corpus, "best", allInCorpus[0][0],\
                #"10p", tenthPerc[-1][0]

                cc = spearmanr([x[0] for x in allInCorpus],
                                [x[1] for x in allInCorpus])

                cc = cc[0]

                cc10 = spearmanr([x[0] for x in tenthPerc],
                              [x[1] for x in tenthPerc])
                cc10 = cc10[0]

                #print cc, cc10
                ccs.append(cc)
                ccs10.append(cc10)

            print stat, "mean", sum(ccs)/len(ccs), "10", sum(ccs10)/len(ccs10)
            print        

    if 1:
        import pylab
        from numpy import array

        stat = "f"

        lst = []
        captions = []
        for solvers,props in evals[4].items():
            captions.append(solvers)

            if len(props) == 1:
                props = 100 * props

            lst.append([prop.stats[stat] for prop in props])

        #print lst

        mat = array(lst)

        import cPickle as pickle
        ff = file('fstats', 'w')
        pickle.dump(captions, ff)
        pickle.dump(mat, ff)

        pylab.boxplot(mat.transpose())
        pylab.xticks(range(1, len(captions)+1), captions)
        pylab.show()

    if 0:
        import pylab
        from numpy import array
        from itertools import cycle

        stat = "121"

        color = cycle(["r","c","y","k","b","g","m"])

        captions = []
        for solvers,props in evals[3].items():
            objs = [prop.stats[solverObj(solvers)] for prop in props]
            stats = [prop.stats[stat] for prop in props]
            captions.append(solvers)
            pylab.scatter(objs, stats, c=color.next())

        pylab.legend(captions)
        pylab.show()
        
