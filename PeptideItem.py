class PeptideItem(object):
    def __init__(self):
        self.mass = 0
        self.peptide = ""
        self.score = 0

    def printPeptide(self):
        print("peptide = " + self.peptide)
        print("mass = " + str(self.mass))
        print("score = " + str(self.score))
