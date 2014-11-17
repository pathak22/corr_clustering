from __future__ import division
import optparse
import sys

from math import log, exp
#from numpy.random.mtrand import dirichlet, multinomial

def addOptions(parseopts):
    parseopts.add_option("-n", "--number", dest="num", default=1,
                         type="int",
                         help="number of points to generate")

if __name__ == "__main__":
    parseopts = optparse.OptionParser()
    addOptions(parseopts)

    (options, args) = parseopts.parse_args()
    num = options.num


    print 1 #dummy training data
    
    print #blank line

    for cluster in range(num):
        print cluster
