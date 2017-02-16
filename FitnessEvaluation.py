from copy import copy,deepcopy
from PeptideItem import PeptideItem
from Score import Score
from Utility import Utility

class FitnessEvaluation(object):

    def __init__(self,method,logIntervalIn,parentMassIn, scoreObjIn, localSearchValidRangeIn):
        self.logInterval = logIntervalIn
        self.fitnessEvaluationCount = 0
        self.scoreStr = ""
        self.fitnessStr = ""
        self.method = method
        self.nextLogCounter = logIntervalIn
        self.maxScoredPeptide = PeptideItem()
        self.utilityObj = Utility()
        self.parentMass = parentMassIn
        self.scoreObj = scoreObjIn
        self.localSearchValidRange = localSearchValidRangeIn

    def increaseFE(self, newPeptide):
        self.fitnessEvaluationCount += 1
        if newPeptide.mass == 0:
            newPeptide.mass = self.utilityObj.calculateMassOfPeptide(newPeptide.peptide)
        if self.method == 3:
            if newPeptide.mass in range(self.parentMass-self.localSearchValidRange,
            self.parentMass+self.localSearchValidRange):
                score = self.scoreObj.sequenceScoreExt(newPeptide.peptide)
                if score > self.maxScoredPeptide.score:
                    self.maxScoredPeptide = copy(newPeptide)
                    self.maxScoredPeptide.score = score
        else:
            if newPeptide.mass == self.parentMass:
                score = self.scoreObj.sequenceScoreExt(newPeptide.peptide)
                if score > self.maxScoredPeptide.score:
                    self.maxScoredPeptide = copy(newPeptide)
                    self.maxScoredPeptide.score = score

        if self.fitnessEvaluationCount >= self.nextLogCounter :
            self.scoreStr = self.scoreStr + str(self.maxScoredPeptide.score) + ";"
            self.fitnessStr = self.fitnessStr + str(self.fitnessEvaluationCount) + ";"
            self.nextLogCounter += self.logInterval


    def getFitnessEvaluationCount(self):
        return self.fitnessEvaluationCount

    def getMaxScore(self):
        return self.maxScoredPeptide.score

    def printFE(self):
        print("FitnessEvaluation count  for method = " + str(self.method))
        print(self.scoreStr)
        print(self.fitnessStr)
