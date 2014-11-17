from __future__ import division

class Scores:
    def __init__(self):
        self.truePos = 0
        self.trueNeg = 0
        self.totalPos = 0
        self.totalNeg = 0
        self.proposed = 0

    def add(self, gold, prop):
        assert(gold in [0,1] and prop in [0,1])

        if gold == 0:
            self.totalNeg += 1

            if prop == gold:
                self.trueNeg += 1
        else:
            self.totalPos += 1

            if prop == gold:
                self.truePos += 1

        if prop == 1:
            self.proposed += 1

    def __repr__(self):
        return "P: %2.2f R: %2.2f F: %2.2f Acc: %2.2g (%d/%d)" % \
               ( 100 * self.prec(), 100 * self.rec(), 100 * self.f(),
                 100 * self.acc(),
                 (self.truePos + self.trueNeg),
                 (self.totalPos + self.totalNeg))

    def __add__(self, other):
        res = Scores()
        res.truePos = self.truePos + other.truePos
        res.trueNeg = self.trueNeg + other.trueNeg
        res.totalPos = self.totalPos + other.totalPos
        res.totalNeg = self.totalNeg + other.totalNeg
        res.proposed = self.proposed + other.proposed

        return res

    def prec(self):
        try:
            return self.truePos / self.proposed
        except ZeroDivisionError:
            return 0

    def rec(self):
        try:
            return self.truePos / self.totalPos
        except ZeroDivisionError:
            return 0

    def acc(self):
        try:
            return ((self.truePos + self.trueNeg) /
                    (self.totalPos + self.totalNeg))
        except ZeroDivisionError:
            return 0

    def f(self):
        den = (self.prec() + self.rec())
        if den == 0:
            return 0
        return (2 * self.prec() * self.rec()) / den        
