from path import path
import sys
import re

import cPickle as pickle

from AIMA import DefaultDict

from nltk_lite.tokenize import regexp, WORDPUNCT
from nltk_lite.stem import Porter

STEMMER = Porter()
def stem(word):
    return STEMMER.stem(word)

class NewsgroupPost:
    def __init__(self, group, vocab):
        self.counts = DefaultDict(0)
        self.subjCounts = DefaultDict(0)
        self.vocab = vocab
        self.group = group

    def read(self, ff):
        for line in ff:
            if re.match("^\S+:", line):
                if line.startswith("Subject:"):
                    line = line.strip("Subject:")
                    self.addCounts(line, subject=True)
                continue

            self.addCounts(line)

    def addCounts(self, line, subject=False):
        #stem(x.lower())
        #words = [x for x in regexp(line, WORDPUNCT)]
        words = [x for x in regexp(line, pattern=r'\w+', gaps=False)]

        for word in words:
            vci = self.vocab.setdefault(word, len(self.vocab))

            self.counts[vci] += 1
            if subject:
                self.subjCounts[vci] += 1

def readNewsgroups(ngDir):
    vocab = {}

    groups = {}

    for group in ngDir.dirs():
        gName = group.basename()
        print >>sys.stderr, "reading", gName

        groups[gName] = []

        for ff in group.files():
            np = NewsgroupPost(gName, vocab)
            np.read(file(ff))

            groups[gName].append(np)

    return vocab,groups

if __name__ == "__main__":
    vocab,groups = readNewsgroups(
        path("DIRECTORY WITH MINI_NEWSGROUPS"))

    print >>sys.stderr, "counting term/doc frequencies"

    wfreqs = DefaultDict(0)
    dfs = DefaultDict(0)

    for group,posts in groups.items():
        for post in posts:
            for word,ct in post.counts.items():
                wfreqs[word] += ct
                dfs[word] += 1

    print >>sys.stderr, "dumping"

    output = path("FILE WHERE WE DUMP THE WORD FREQS")
    fh = file(output, 'w')
    pickle.dump(vocab, fh)
    pickle.dump(groups, fh)
    pickle.dump(wfreqs, fh)
    pickle.dump(dfs, fh)
