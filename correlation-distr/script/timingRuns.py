from Hogwash import Session
from Hogwash.Results import ResultsFile
from waterworks.Processes import bettersystem

import sys
from path import path
import os
from StringIO import StringIO

def on64():
    arch = os.uname()[-1]
    return arch == "x86_64"

def myBin(mydir):
    upOne = mydir.parent

    if on64() and False: #no 64-bit support right now!
        res = upOne/"bin64"
    else:
        res = upOne/"bin32"
    print res
    assert(res.exists())
    return res

def runtask(models, datadir):
    mfile = datadir/"matrix"
    bin = myBin(datadir)
    solver = bin/"chainedSolvers"
    mfile = datadir/"matrix"
    outfile = datadir/("output-%s" % "-".join(models[0].split()))

    output = StringIO()
    cmd = "%s %s %s" % (solver, " ".join(models), mfile)
    print "Running", cmd
    bettersystem(cmd, stdout=output)
    print "Done, postprocessing output..."
    outfh = file(outfile, 'w')
    for line in output.getvalue().split("\n"):
        if not line.strip().isdigit():
            print line
        else:
            print  >>outfh, line
    outfh.close()
    
    return ResultsFile(outfile)

if __name__ == "__main__":
    workdir = path(sys.argv[1])

    print "Working directory", workdir
    if not workdir.exists():
        workdir.mkdir()

        print "No data... generating data now"

        for power in range(1,4):
            for layer in range(9):
                dsize = (layer + 1) * 10**power
                print "Data: size", dsize
                if dsize >= 5000:
                    break

                for run in range(3):
                    print "Run:", run

                    ddir = workdir/("data%d-%d" % (dsize, run))
                    ddir.mkdir()

                    output = StringIO()
                    cmd = "python script/genGaussians.py -n %d -t 10 -f 3 -a 1"
                    bettersystem(cmd % dsize, stdout=output)

#                    bettersystem("python script/genFeatureless.py -n %d" %
#                                 dsize, stdout=output)

                    dfileName = ddir/"data"
                    dfile = file(dfileName, 'w')
                    dfile.write(output.getvalue()),
                    dfile.close()

                    output = StringIO()
                    bettersystem("python script/classify.py -c max-ent %s" %
                                 dfileName, stdout=output)

#                    bettersystem("python script/classify.py %s" %
#                                 dfileName, stdout=output)
                    mfileName = ddir/"matrix"
                    mfile = file(mfileName, 'w')
                    mfile.write(output.getvalue()),
                    mfile.close()

    models = [
        ["first"],
        ["best"],
        ["vote"],
        ["pivot"],
        ["boem"],
        ["first boem"],
        ["best boem"],
        ["vote boem"],
        ["pivot boem"],
        ["lazylp first"],
        ["lazylp best"],
        ["lazylp vote"],
        ["lazylp pivot"],
        ["lazylp first boem"],
        ["lazylp best boem"],
        ["lazylp vote boem"],
        ["lazylp pivot boem"],
        ["ilp"],
        ["sumsdp vote"],
        ]

    args = []
    for model in models:
        for datadir in workdir.dirs():
            if "bin" in datadir.basename():
                continue
            args.append((model, datadir.abspath()))

    ssn = Session("timingRuns", "runtask", args,
                  name=("hog%s" % workdir.basename()))
