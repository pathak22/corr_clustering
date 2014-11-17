from __future__ import division

import optparse
import sys

from iterextras import batch
from itertools import izip, repeat, count

from scores import Scores

from Classifiers import Point, readDataFile
from ExtraClassifiers import RandomErrClassifier, NaiveBayesClassifier
from maxEntDots import MaxEntDots
from maxEntNewsgroups import MaxEntNewsgroups

from StringIO import StringIO

from Hogwash import Session
from Hogwash.Sty import Sty
from time import sleep

def partLabel(classifier, rows):
    out = StringIO()
    classifier.partialLabel(rows, out)
    return (out.getvalue(), classifier.scores)

if __name__ == "__main__":
    parseopts = optparse.OptionParser()
    parseopts.add_option("-c", "--classifier", dest="classifier",
                         default="random-errors",
                         help="which classifier to use")
    parseopts.add_option("-a", "--all-classifiers", dest="list", default=False,
                         action="store_true",
                         help="list classifiers available")
    parseopts.add_option("-s", "--simple", dest="simple", default=False,
                         action="store_true",
                         help="effectively sets aIn=bOut=1, aOut=bIn=0 for random errors")
    parseopts.add_option("", "--aIn", dest="aIn", default=4,
                         type="float",
                         help="aIn parameter for random errors")
    parseopts.add_option("", "--bIn", dest="bIn", default=1,
                         type="float",
                         help="bIn parameter for random errors")
    parseopts.add_option("", "--aOut", dest="aOut", default=1,
                         type="float",
                         help="aOut parameter for random errors")
    parseopts.add_option("", "--bOut", dest="bOut", default=4,
                         type="float",
                         help="bOut parameter for random errors")
    parseopts.add_option("-e", "--epsilon", dest="epsilon", default=.1,
                         type="float",
                         help="""random errors are made on an epsilon fraction
                         of edges""")
    parseopts.add_option("-n", "--newsgroups", dest="newsgroups", default=
                         "FILE WHERE WE DUMP THE WORD FREQS",
                         type="string",
                         help="""dump file with newsgroups data stats""")
    parseopts.add_option("", "--check", dest="check", default=False,
                         action="store_true",
                         help="""do a randomized check of the classifier
                         instead of generating a matrix""")
    parseopts.add_option("-d", "--distributed", dest="distributed",
                         default="",
                         type="string",
                         help="""use hogwash to run the labeling step
                         in parallel; argument is the name of the session""")

    (options, args) = parseopts.parse_args()

    CLASSIFIERS = {
        "random-errors" : RandomErrClassifier,
        "naive-bayes" : NaiveBayesClassifier,
        "max-ent" : MaxEntDots,
        "mxnews" : MaxEntNewsgroups,
        }

    if options.list:
        for clType,cl in CLASSIFIERS.items():
            print clType, "\n", cl.__doc__
        sys.exit(0)

    dataFile = file(args[0])
    print >>sys.stderr, "reading", dataFile

    trainPoints = readDataFile(dataFile)
    print >>sys.stderr, "read", len(trainPoints), "training points"

    testPoints = readDataFile(dataFile)
    print >>sys.stderr, "read", len(testPoints), "testing points"

    try:
        classifier = CLASSIFIERS[options.classifier](
            trainPoints, testPoints, options)
    except KeyError:
        print "No such classifier as", options.classifier
        sys.exit(1)

    print >>sys.stderr, "training..."
    classifier.train()

    if options.check:
        print >>sys.stderr, "evaluating..."
        classifier.check()
    elif options.distributed:
        print >>sys.stderr, "distributed labeling..."
        args = list(izip(repeat(classifier),
                         batch(range(len(testPoints)), 30)))
        ssn = Session("classify", "partLabel", args, name=options.distributed)
        print >>sys.stderr, "session created..."

        if 0:
            sty = Sty("workingSession")
            sty.do_clear("all")
            sty.do_runlocally("")
            jobs_done = False
            while not jobs_done:
                stats = ssn.get_job_status_report()
                queued = []
                if "queued" in stats:
                    queued = stats["queued"]
                started = []
                if "started" in stats:
                    started = stats["started"]
                if not (queued + started):
                    jobs_done = True
                else:
                    print >>sys.stderr,\
                          len(started), "started", len(queued), "queued"
                    sleep(10)
                if "error" in stats and stats["error"]:
                    print >>sys.stderr, "Errors occurred", stats["error"]
                    sys.exit(1)

            finalScores = Scores()
            print len(testPoints)
            for job in ssn:
                res, scores = job.results
                print res,
                finalScores += scores
            classifier.scores = finalScores
        else:
            sys.exit(0)

    else:
        print >>sys.stderr, "labeling..."
        classifier.label()

    print >>sys.stderr, "test classifier performance"
    print >>sys.stderr, classifier.scores
