import cPickle as pickle
import sys
import pylab
from numpy import array, concatenate

if __name__ == "__main__":
    stat = "f"
    ff1 = file('%sstats' % stat)
    addSolvers = pickle.load(ff1)
    addMat = pickle.load(ff1)

    ff2 = file('log%sstats' % stat)
    logSolvers = pickle.load(ff2)
    logMat = pickle.load(ff2)

    print addSolvers, logSolvers

    solvers = addSolvers + logSolvers
    print solvers
    mat = concatenate((addMat, logMat))
    print mat.shape

    solvers = array([" ".join(ss) for ss in solvers])

    #let's snip out some solvers that suck
    SUCK = ["first", "log first", "best", "log best", "pivot",
            "log pivot", "log sa"]
    keep = []
    for ct,solver in enumerate(solvers):
        if solver not in SUCK:
            keep.append(ct)

    newmat = mat[keep, :]
    solvers = solvers[keep]

    solvers = [solver.replace("log ", "L").replace(" boem", "B") for
               solver in solvers]

    pylab.boxplot(newmat.transpose(), )
    pylab.xticks(range(1, len(solvers)+1), solvers, rotation=-45)
    pylab.show()

