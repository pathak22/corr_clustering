from Hogwash.Results import ResultsFile #type for file created by hw job
from Hogwash.Action import Action #supertype for runnable objects
from Hogwash.Errors import BadExitCode #error if the program crashed
#get a job-specific name
from Hogwash.Helpers import make_job_output_filename, get_cpu_bitness
from waterworks.Processes import bettersystem #run an external command
from evaluate import readClMat, Eval
from Classifiers import readDataFile

from analysis.chatStats import readChat, copyChat, nonSys
from analysis.distMat import getLoc, get121, getF

from path import path
from StringIO import StringIO #store output of process
from itertools import izip

#class representing an experiment to run
class Experiment(Action):
    def __init__(self, solverCmd, solvers, mfile):
        self.solverCmd = solverCmd.abspath()
        self.solvers = solvers
        self.mfile = mfile.abspath()

    def run(self, hogwash_job):
        on64 = (get_cpu_bitness() == 64)
        if on64:
            par = self.solverCmd.parent
            self.solverCmd = par/"../bin64"/self.solverCmd.basename()

        outfile = make_job_output_filename(hogwash_job, "output")
        cmd = "%s %s %s > %s" % (self.solverCmd, " ".join(self.solvers),
                                 self.mfile, outfile)

        print cmd

        status = bettersystem(cmd)
        if status != 0:
            raise BadExitCode(status)

        return ResultsFile(outfile)

class Evaluation(Action):
    def __init__(self, exp, gold, mat):
        self.exp = exp
        self.gold = gold.abspath()
        self.mat = mat.abspath()

    def run(self, hogwash_job):
        self.replace_actions(hogwash_job) #bug in hogwash
        proposal = readDataFile(file(self.exp))
        goldFile = file(self.gold)
        readDataFile(goldFile) #discard training data
        gold = readDataFile(goldFile)
        matrix = readClMat(file(self.mat))
        evaluator = Eval(gold, proposal, matrix, self.exp)
        return evaluator

class ChatEvaluation(Evaluation):
    def __init__(self, exp, nodeVectorGold, chatGolds, mat):
        Evaluation.__init__(self, exp, nodeVectorGold, mat)
        self.chatGolds = [x.abspath() for x in chatGolds]

    def run(self, hogwash_job):
        res = Evaluation.run(self, hogwash_job)

        golds = []
        for chat in self.chatGolds:
            print "Reading", chat
            golds.append([x for x in readChat(chat)])

        relabel = copyChat(golds[0])
        proposal = readDataFile(file(self.exp)) #opens twice, whatever
        assert(len(proposal) == len(nonSys(relabel)))
        for line,prop in izip(nonSys(relabel), proposal):
            line.thread = prop

        scores = {}
        for metric in [get121, getF, getLoc]:
            print "Computing", metric.__name__
            mscores = [metric(gold, relabel) for gold in golds]
            scores[metric.__name__] = mscores

        res.chatScores = scores
        return res
