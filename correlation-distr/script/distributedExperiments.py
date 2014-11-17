from Hogwash import Session #main hogwash class
from Hogwash.Results import ResultsFile #type for file created by hw job
from Hogwash.Action import Action #supertype for runnable objects
from Hogwash.Errors import BadExitCode #error if the program crashed
from waterworks.Processes import bettersystem #run an external command

import sys
from path import path
import os
from shutil import copy
from iterextras import batch

from StringIO import StringIO #store output of process

from DistributedExperiment import Experiment, Evaluation

if __name__ == "__main__":
    sessionName = sys.argv[1]
    matrices = [path(x) for x in sys.argv[2:]]

    sessionBinDir = path(sessionName)/"bin32"
    solverCmd = sessionBinDir/"chainedSolvers"

    args = []

    greedyMethods = [["first"], ["best"], ["vote"], ["pivot"], ["sa"]]
    identMethods = [["%s-id" % x[0]] for x in greedyMethods]
    greedyMethods = 100 * greedyMethods
    greedyBOEM = [ge+["boem"] for ge in greedyMethods]

    identBOEM =  [ge+["boem"] for ge in identMethods]
    identMethods += identBOEM

    for solvers in greedyMethods + greedyBOEM + [["sdp2l", "pivot", "boem"],
                                                 ["boem"]] + identMethods:
        for matrix in matrices:
            #exp = Experiment(solverCmd, ["log"] + solvers, matrix)
            exp = Experiment(solverCmd, solvers, matrix)
            args.append(exp)
            ev = Evaluation(exp, matrix.stripext()+".truth", matrix)
            args.append(ev)

    session = Session("Hogwash.Action", "action_runner", args, name=sessionName)

    deps = {}
    for expt,evalt in batch(session, 2):
            deps[evalt] = expt
    session.set_dependencies(deps)

    if not sessionBinDir.exists():
        sessionBinDir.mkdir()

    copy("bin32/chainedSolvers", sessionBinDir)

    session64BinDir = path(sessionName)/"bin64"
    if not session64BinDir.exists():
        session64BinDir.mkdir()

    copy("bin64/chainedSolvers", session64BinDir)
