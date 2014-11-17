from megam import Megam
from Classifiers import Classifier

class MaxEntDots(Classifier):
    """This classifier runs logistic regression on the absolute coordinate
    differences between the input points."""
    def __init__(self, train, label, extra):
        Classifier.__init__(self, train, label, extra)
        self.classifier = Megam()

    def makeFeatDict(self, p1, p2):
        fdict = {}
        for feat in set(p1.feats.keys() + p2.feats.keys()):
            val1 = p1.feats.get(feat, 0)
            val2 = p2.feats.get(feat, 0)

            diff = abs(val1 - val2)
            fdict["DIFF_%g" % feat] = diff
        return fdict

    def addTrainingInst(self, instance):
        (label, p1, p2) = instance
        fdict = self.makeFeatDict(p1, p2)
        self.classifier.add(label, fdict)

    def estimate(self):
        self.classifier.train()

    def classify(self, instance):
        (label, p1, p2) = instance
        fdict = self.makeFeatDict(p1, p2)

        res = self.classifier.classify(fdict)

        self.scores.add(label, (res > .5))
        return res
