from copy import copy,deepcopy
import numpy as np
from FitnessEvaluation import FitnessEvaluation
from Score import Score
from Init import Init
from Operator import Operator
from PeptideItem import PeptideItem
from Utility import Utility

class PeptideSequence(object):
    def __init__(self, initObjIn, methodIn, acoIterationIn, finalLSIterationIn, AcoLocalFitnessEvaluation,
     LogInterval):
        self.AcoLocalFitnessEvaluation = AcoLocalFitnessEvaluation
        self.initObj = initObjIn
        self.localSearchValidRange = 10
        self.method = methodIn
        self.scoreObj = Score(initObjIn.spectrumMap)
        self.fitnessObj = FitnessEvaluation(methodIn, LogInterval, initObjIn.parentMass, self.scoreObj,
        self.localSearchValidRange)
        self.operatorObj = Operator(self.fitnessObj, self.scoreObj)
        self.pheromone = initializePheromone(initObjIn.aminoAcidMass, initObjIn.parentMass)
        self.utilityObj = Utility()
        self.highestScoredPeptideMap = {}
        self.acoIteration = acoIterationIn
        self.finalLSIteration = finalLSIterationIn
        self.localSearchIteration = 1


    def getAcoLocalFitnessEvaluation(self):
        return self.AcoLocalFitnessEvaluation

    def executeMethod(self):
        if self.method == 1:
            return self.acoLocal()
        elif self.method == 2:
            return self.acoOnly()
        elif self.method == 3:
            return self.LocalSearchOnly()
        elif self.method == 4:
            return self.BnB()



    def acoLocal(self):
        self.AcoLocalFitnessEvaluation = 0
        print("in AcoLocal")
        print("acoIteration = " + str(self.acoIteration))
        print("finalLSIteration = " + str(self.finalLSIteration))
        #aco
        for epoch in range(2):
            self.pheromone.fill(1.)
            probability = self.pheromone.copy()
            probability = probability + 0.1
            probability = probability/probability.sum(axis=1)[:,None]
            for k in range(self.acoIteration):
                peptideItemList = []
                #print("no of iteration = " + str(k))
                for i in range(self.initObj.ants):
                    newPeptide = createAntPath(probability,False,self.initObj, self.scoreObj, self.utilityObj)
                    peptideItemList.append(newPeptide)
                    self.fitnessObj.increaseFE(newPeptide)

                modifiedPeptideItemList = performLocalSearch(peptideItemList,self.localSearchIteration,
                self.operatorObj, self.initObj)
                self.pheromone = calculatePheromone(modifiedPeptideItemList, self.initObj, self.pheromone)
                saveHighestScoredPeptides(modifiedPeptideItemList, False, self.initObj,
                self.highestScoredPeptideMap, self.localSearchValidRange)
                #print("highestScoredPeptideMap within iterations")
                #printPeptideMap(self.highestScoredPeptideMap)
                probability = copy(self.pheromone)
                probability = probability + 0.1
                probability = probability/probability.sum(axis=1)[:,None]

        print("HighestScoredPeptideMap after first " + str(self.acoIteration)  + "  iterations")
        printPeptideMap(self.highestScoredPeptideMap)
        peptideItemList = getHighestScoredPeptides(True, self.highestScoredPeptideMap, self.scoreObj)
        print("highestScored peptides after conversion to SequenceScoreExt")
        for i in range(len(peptideItemList)):
            print(peptideItemList[i].peptide + "-> (" + str(peptideItemList[i].score) + ")" )

        self.highestScoredPeptideMap.clear()
        peptideItemList = localSearchExt(peptideItemList, False, self.finalLSIteration, self.method,
        self.AcoLocalFitnessEvaluation, self.operatorObj, self.initObj, self.highestScoredPeptideMap,
        self.scoreObj, self.fitnessObj, 0)
        self.AcoLocalFitnessEvaluation = self.fitnessObj.getFitnessEvaluationCount()
        printEvaluationCounters(self.AcoLocalFitnessEvaluation, self.fitnessObj.getFitnessEvaluationCount())
        self.fitnessObj.printFE()
        print("HighestScoredPeptideMap after localSearchExt "+ str(self.finalLSIteration) + " iterations")
        printPeptideMap(self.highestScoredPeptideMap)
        return self.highestScoredPeptideMap


    def acoOnly(self):
        print("AcoOnly")
        print("acoIteration = " + str(self.acoIteration))
        print("finalLSIteration = " + str(self.finalLSIteration))
        #aco
        oneFourthMileStone = int(self.AcoLocalFitnessEvaluation/4)
        halfMileStone = int(self.AcoLocalFitnessEvaluation / 2)
        threeFourthMileStone = int(self.AcoLocalFitnessEvaluation * 3/4)
        restart = False
        reachedOneFourth = False
        reachedThreeFourth = False
        peptideMapFirstEpoch = {}

        for epoch in range(2):
            #print("epoch = " + str(epoch))
            self.pheromone.fill(1.)
            probability = self.pheromone.copy()
            probability = probability + 0.1
            probability = probability/probability.sum(axis=1)[:,None]
            extensiveScoring = False

            while self.fitnessObj.getFitnessEvaluationCount() < self.AcoLocalFitnessEvaluation:
                peptideItemList = []
                #print("no of iteration = " + str(k))
                for i in range(self.initObj.ants):
                    newPeptide = createAntPath(probability, extensiveScoring, self.initObj, self.scoreObj, self.utilityObj)
                    peptideItemList.append(newPeptide)
                    self.fitnessObj.increaseFE(newPeptide)

                self.pheromone = calculatePheromone(peptideItemList, self.initObj, self.pheromone)
                saveHighestScoredPeptides(peptideItemList, False, self.initObj,
                self.highestScoredPeptideMap, self.localSearchValidRange)
                probability = self.pheromone.copy()
                probability = probability + 0.1
                probability = probability/probability.sum(axis=1)[:,None]

                if reachedOneFourth == False:
                    if self.fitnessObj.getFitnessEvaluationCount() >= oneFourthMileStone:
                        extensiveScoring = True
                        reachedOneFourth = True
                        peptideItemList = getHighestScoredPeptides(True, self.highestScoredPeptideMap,
                        self.scoreObj)
                        self.highestScoredPeptideMap.clear()
                        saveHighestScoredPeptides(peptideItemList, False, self.initObj,
                        self.highestScoredPeptideMap, self.localSearchValidRange)

                elif restart == False:
                    if self.fitnessObj.getFitnessEvaluationCount() >= halfMileStone:
                        restart = True
                        peptideMapFirstEpoch = deepcopy(self.highestScoredPeptideMap)
                        self.highestScoredPeptideMap.clear()
                        break
                elif reachedThreeFourth == False:
                    if self.fitnessObj.getFitnessEvaluationCount() >= threeFourthMileStone:
                        extensiveScoring = True
                        reachedThreeFourth = True
                        peptideItemList = getHighestScoredPeptides(True, self.highestScoredPeptideMap,
                        self.scoreObj)
                        self.highestScoredPeptideMap.clear()
                        saveHighestScoredPeptides(peptideItemList, False, self.initObj,
                        self.highestScoredPeptideMap, self.localSearchValidRange)

        printEvaluationCounters(self.AcoLocalFitnessEvaluation, self.fitnessObj.getFitnessEvaluationCount())
        self.fitnessObj.printFE()
        peptideItemList = getHighestScoredPeptides(True, peptideMapFirstEpoch,self.scoreObj)
        saveHighestScoredPeptides(peptideItemList, False, self.initObj,
        self.highestScoredPeptideMap, self.localSearchValidRange)

        print("highestScoredPeptideMap in acoOnly after"+ str(self.fitnessObj.getFitnessEvaluationCount()) + " fitness evaluation")
        printPeptideMap(self.highestScoredPeptideMap)
        return self.highestScoredPeptideMap


    def LocalSearchOnly(self):
        print("LocalSearchOnly")
        print("acoIteration = " + str(self.acoIteration))
        print("finalLSIteration = " + str(self.finalLSIteration))

        randomlyGeneratedPeptide = []
        for i in range(self.initObj.ants):
            generatedPeptide = PeptideItem()
            while generatedPeptide.mass <= self.initObj.parentMass:
                rand = np.random.randint(0,len(self.initObj.aminoAcidMass))
                generatedPeptide.peptide += str(self.initObj.aminoAcidMass[rand]) + "-"
                generatedPeptide.mass += self.initObj.aminoAcidMass[rand]
            generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
            generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)
            randomlyGeneratedPeptide.append(generatedPeptide)
            self.fitnessObj.increaseFE(generatedPeptide)


        peptideItemList = copy(randomlyGeneratedPeptide)
        for counter in range(2):
            for k in range(self.acoIteration):
                if len(peptideItemList) > 0:
                    modifiedPeptideItemList = performLocalSearch(peptideItemList,self.localSearchIteration,
                    self.operatorObj, self.initObj)
                    saveHighestScoredPeptides(modifiedPeptideItemList, True,self.initObj,
                    self.highestScoredPeptideMap, self.localSearchValidRange)
                    peptideItemList = getHighestScoredPeptides(False, self.highestScoredPeptideMap,
                    self.scoreObj)
                else:
                    peptideItemList = copy(randomlyGeneratedPeptide)
                    break


        peptideItemList = getHighestScoredPeptides(True, self.highestScoredPeptideMap, self.scoreObj)
        self.highestScoredPeptideMap.clear()

        peptideItemList = localSearchExt(peptideItemList, True, self.finalLSIteration, self.method,
        self.AcoLocalFitnessEvaluation, self.operatorObj, self.initObj, self.highestScoredPeptideMap,
        self.scoreObj, self.fitnessObj, self.localSearchValidRange)

        printEvaluationCounters(self.AcoLocalFitnessEvaluation, self.fitnessObj.getFitnessEvaluationCount())
        self.fitnessObj.printFE()
        print("highestScoredPeptideMap in LocalSearchOnly after " +
        str(self.fitnessObj.getFitnessEvaluationCount()) + " fitness evaluation")
        printPeptideMap(self.highestScoredPeptideMap)
        return self.highestScoredPeptideMap


    def BnB(self):
        leaderboard = []
        leaderboard.append("")
        leaderboard = expand(leaderboard, self.initObj.aminoAcidMass)
        leaderPeptide = ""
        maxScore = 0
        counter = 0
        maxScoredPeptide = []
        generatedPeptide = PeptideItem()

        while (self.fitnessObj.getFitnessEvaluationCount()<self.AcoLocalFitnessEvaluation) and (len(leaderboard)>0):
            leaderboard = expand(leaderboard, self.initObj.aminoAcidMass)
            tmpleaderboard =  deepcopy(leaderboard)
            for i in range(len(leaderboard)):
                peptide = leaderboard[i]
                mass = self.utilityObj.calculateMassOfPeptide(peptide)
                if mass == self.initObj.parentMass:
                    cyscore = self.scoreObj.cycloScore(peptide)
                    leaderPeptideCyScore = self.scoreObj.cycloScore(leaderPeptide)
                    generatedPeptide.peptide = peptide
                    generatedPeptide.score = cyscore
                    self.fitnessObj.increaseFE(generatedPeptide)
                    if leaderPeptide != "":
                        generatedPeptide.peptide = leaderPeptide
                        generatedPeptide.score = leaderPeptideCyScore
                        self.fitnessObj.increaseFE(generatedPeptide)

                    if cyscore >= leaderPeptideCyScore:
                        leaderPeptide = peptide
                        if cyscore > maxScore:
                            maxScoredPeptide = []
                            maxScoredPeptide.append(peptide)
                        elif cyscore == maxScore:
                            maxScoredPeptide.append(peptide)
                        maxScore = cyscore
                elif mass > self.initObj.parentMass:
                    tmpleaderboard.remove(peptide)

            if len(tmpleaderboard) > 0:
                leaderboard = trim(tmpleaderboard, self.initObj.n, self.fitnessObj, self.scoreObj)

            else:
                leaderboard = tmpleaderboard
            counter += 1
        print("counter = " + str(counter))

        modifiedPeptideItemList = []
        for i in range(len(maxScoredPeptide)):
            score = self.scoreObj.sequenceScoreExt(maxScoredPeptide[i])
            generatedPeptide = PeptideItem()
            generatedPeptide.peptide = maxScoredPeptide[i]
            generatedPeptide.score = score
            generatedPeptide.mass = self.utilityObj.calculateMassOfPeptide(maxScoredPeptide[i])
            modifiedPeptideItemList.append(generatedPeptide)
        saveHighestScoredPeptides(modifiedPeptideItemList, False, self.initObj,
        self.highestScoredPeptideMap, 0)

        printEvaluationCounters(self.AcoLocalFitnessEvaluation, self.fitnessObj.getFitnessEvaluationCount())
        self.fitnessObj.printFE()

        print("HighestScoredPeptideMap in BnB")
        printPeptideMap(self.highestScoredPeptideMap)
        return self.highestScoredPeptideMap


def expand(peptides, aminoAcidMass):
    extendedPeptide = []
    if peptides[0] == "":
        extendedPeptide = aminoAcidMass
    else:
        for i in range(len(peptides)):
            for j in range(len(aminoAcidMass)):
                tmpPeptide = str(peptides[i]) + '-' +str(aminoAcidMass[j])
                extendedPeptide.append(tmpPeptide)
    return extendedPeptide


def trim(leaderboard, n, fitnessObj, scoreObj):
    leaderboardLen = len(leaderboard)
    linearScoreMap = {}
    linearScoreList = []
    generatedPeptide = PeptideItem()
    for i in range(leaderboardLen):
        peptide = leaderboard[i]
        score = scoreObj.linearScore(peptide)
        generatedPeptide.peptide = peptide
        generatedPeptide.score = score
        fitnessObj.increaseFE(generatedPeptide)
        linearScoreList.append(score)
        if score not in linearScoreMap:
            linearScoreMap[score] = []
        linearScoreMap[score].append(peptide)
    linearScoreList.sort(reverse = True)
    for j in range(n,len(linearScoreList)):
        if linearScoreList[j] < linearScoreList[n-1]:
            del linearScoreList[(j):]
            break

    leaderboard = []
    peptides = linearScoreMap[linearScoreList[0]]
    for j in range(len(peptides)):
        leaderboard.append(peptides[j])
    for i in range(1,len(linearScoreList)):
        if linearScoreList[i] == linearScoreList[i-1]:
            continue
        peptides = linearScoreMap[linearScoreList[i]]
        for j in range(len(peptides)):
            leaderboard.append(peptides[j])
    return leaderboard



def localSearchExt(peptideItemList, isLocalSearchOnly, finalLSIteration, method, AcoLocalFitnessEvaluation,
operatorObj, initObj, highestScoredPeptideMap, scoreObj, fitnessObj, localSearchValidRange):
    modifiedPeptideItemList = copy(peptideItemList)
    j = 0
    counter = finalLSIteration
    if method > 1:
        j = fitnessObj.getFitnessEvaluationCount()
        counter = AcoLocalFitnessEvaluation
    #print("j = " + str(j))
    #print("AcoLocalFitnessEvaluation = " + str(AcoLocalFitnessEvaluation))
    #print("counter = " + str(counter))
    #print("FitnessEvaluation = " + str(FitnessEvaluation))
    while j < counter and len(modifiedPeptideItemList) > 0:
        i = 0
        length = len(modifiedPeptideItemList)
        while i < length:
            tmpPeptide = modifiedPeptideItemList[i]
            resultlist = []
            #coin toss for Transposition, swap, fission, fussion will go here
            #1 for mutate
            #2 for Transposition
            #3 for swap
            #4 for Fusion
            #5 for Fission
            #6 for MergeAndBreak
            randomOperator = np.random.choice([1,2,3,4,5,6], p=[0.1,0.3,0.3,0.05,0.05,0.2])
            #print("randomOperator = " + str(randomOperator))
            if randomOperator == 1:
                resultlist = operatorObj.mutate(tmpPeptide,True, initObj)
            elif randomOperator == 2:
                resultlist = operatorObj.transposition(tmpPeptide,True)
            elif randomOperator == 3:
                resultlist = operatorObj.swap(tmpPeptide,True)
            elif randomOperator == 4:
                resultlist = operatorObj.fusion(tmpPeptide,True, initObj.mapToConcatenate)
            elif randomOperator == 5:
                resultlist = operatorObj.fission(tmpPeptide,True, initObj.mapToSplit)
            else:
                resultlist = operatorObj.mergeAndBreak(tmpPeptide,True, initObj.mergeAndBreakMap)
            if len(resultlist) > 0:
                modifiedPeptideItemList.extend(resultlist)
            i += 1

        saveHighestScoredPeptides(modifiedPeptideItemList, isLocalSearchOnly, initObj,
        highestScoredPeptideMap, localSearchValidRange)
        modifiedPeptideItemList = getHighestScoredPeptides(False, highestScoredPeptideMap, scoreObj)
        '''print("modifiedPeptideItemList after path process j = " + str(j))
        for k in modifiedPeptideItemList:
            print(k.peptide + "(" + str(k.score) + " -> "+str(k.mass)+")")'''
        if method > 1:
            j = fitnessObj.getFitnessEvaluationCount()
        else:
            j += 1
    reversePeptideList = operatorObj.reverse(modifiedPeptideItemList,True)
    modifiedPeptideItemList.extend(reversePeptideList)
    saveHighestScoredPeptides(modifiedPeptideItemList, isLocalSearchOnly, initObj,
    highestScoredPeptideMap,localSearchValidRange)
    return modifiedPeptideItemList



def getHighestScoredPeptides(withNewScore, highestScoredPeptideMap, scoreObj):
    highestScoredPeptideList = []
    keys = highestScoredPeptideMap.keys()
    for i in keys:
        values = highestScoredPeptideMap[i]
        for peptide in values:
            if withNewScore:
                peptide.score = scoreObj.sequenceScoreExt(peptide.peptide)
            highestScoredPeptideList.append(peptide)

    return highestScoredPeptideList



def saveHighestScoredPeptides(peptideItemList,isLocalSearchOnly,initObj, highestScoredPeptideMap, localSearchValidRange):
    peptideScoreMap = {}
    for i in range(len(peptideItemList)):
        tmpPeptide = peptideItemList[i]
        if isLocalSearchOnly:
            if tmpPeptide.mass in range(initObj.parentMass-localSearchValidRange, initObj.parentMass+localSearchValidRange):
                if tmpPeptide.score not in peptideScoreMap:
                    peptideScoreMap[tmpPeptide.score] = []
                peptideScoreMap[tmpPeptide.score].append(tmpPeptide)
        else:
            if tmpPeptide.mass == initObj.parentMass:
                if tmpPeptide.score not in peptideScoreMap:
                    peptideScoreMap[tmpPeptide.score] = []
                peptideScoreMap[tmpPeptide.score].append(tmpPeptide)

    mergedKeys = list(highestScoredPeptideMap.keys())
    #print("highestScoredPeptideMap.keys()")
    #print(mergedKeys)
    peptideScoreMapKeys = list(peptideScoreMap.keys())
    #print("peptideScoreMap.keys()")
    #print(peptideScoreMapKeys)

    mergedKeys.extend(peptideScoreMapKeys)
    mergedKeys = np.unique(mergedKeys)
    mergedKeys[::-1].sort()
    #print("after merge and sort")
    #print(mergedKeys)
    count = 0
    tmpMap = copy(highestScoredPeptideMap)
    highestScoredPeptideMap.clear()
    for i in range(len(mergedKeys)):
        k = mergedKeys[i]
        v = tmpMap.get(k,[])
        v.extend(peptideScoreMap.get(k,[]))
        other = {k:[]}
        highestScoredPeptideMap.update(other)
        values = []
        for j in range(len(v)):
            if v[j].peptide not in values:
                values.append(v[j].peptide)
                highestScoredPeptideMap[k].append(v[j])
        count += len(values)
        if count >= initObj.n:
            break


def performLocalSearch(peptideItemList, localSearchIteration, operatorObj, initObj):
    #print("incoming peptideItemList in performLocalSearch")
    #for i in peptideItemList:
        #print(i.peptide + "(" + str(i.score) + " -> "+str(i.mass)+")")
    modifiedPeptideItemList = copy(peptideItemList)
    for m in range(localSearchIteration):
        i = 0
        while i < len(modifiedPeptideItemList):
            tmpPeptide = modifiedPeptideItemList[i]
            resultlist = operatorObj.mutate(tmpPeptide, False, initObj)

            if len(resultlist) > 0:
                del modifiedPeptideItemList[i]
                for k in range(len(resultlist)):
                    modifiedPeptideItemList.insert(i+k,resultlist[k])
                i += len(resultlist)
            else:
                i += 1

    #print("outgoing modifiedPeptideItemList")
    #for k in modifiedPeptideItemList:
        #print(k.peptide + "(" + str(k.score) + " -> "+str(k.mass)+")")
    return modifiedPeptideItemList


def createAntPath(probability,extensiveScoring,initObj,scoreObj, utilityObj):
    newPeptide = PeptideItem()
    while newPeptide.mass <= initObj.parentMass:
        if (initObj.parentMass-newPeptide.mass) < initObj.aminoAcidMass[0]:
            break
        splittedNewPeptide = newPeptide.peptide.split("-")
        peptideLen = len(splittedNewPeptide)
        probabilityRow = probability[peptideLen-1].copy()
        dimensionProbability = probability.shape
        massToBeAdded = np.random.choice(initObj.aminoAcidMass, p=probabilityRow)
        newPeptide.peptide += str(massToBeAdded) + "-"
        newPeptide.mass = utilityObj.calculateMassOfPeptide(newPeptide.peptide)
    #remove last -
    newPeptide.peptide = newPeptide.peptide.strip("-")
    if extensiveScoring:
        newPeptide.score = scoreObj.sequenceScoreExt(newPeptide.peptide)
    else:
        newPeptide.score = scoreObj.sequenceScore(newPeptide.peptide)
    return newPeptide


def calculatePheromone(peptideItemList, initObj, pheromone):
    evaporateQuantity = 1 - initObj.pheromoneConstant
    pheromone = pheromone * evaporateQuantity

    # add pheromone for all ants
    for i in range(len(peptideItemList)):
        if peptideItemList[i].score > 0:
            tmpPeptide = peptideItemList[i].peptide
            aminoAcidsInTmpPeptide = tmpPeptide.split('-')
            aminoAcidsInTmpPeptide = [int(j) for j in aminoAcidsInTmpPeptide]
            pheromoneToBeAdded = peptideItemList[i].score * initObj.q
            for j in range(len(aminoAcidsInTmpPeptide)):
                columnIndexPheromonMatrix = initObj.aminoAcidMass.index(aminoAcidsInTmpPeptide[j])
                rowIndexInPheromoneMatrix = j
                pheromone[rowIndexInPheromoneMatrix,columnIndexPheromonMatrix] += pheromoneToBeAdded
    return pheromone


def initializePheromone(aminoAcidMassIn,parentMassIn):
    noOfColumns = len(aminoAcidMassIn)
    noOfRows = np.ceil(parentMassIn / aminoAcidMassIn[0])
    #print("In Pheromone matrix, number of columns = " + str(noOfColumns) + " no of rows = " + str(noOfRows))
    pheromone = np.ones((noOfRows,noOfColumns),dtype = float)
    pheromoneDimension = pheromone.shape
    #self.rowsCount = pheromoneDimension[0]
    #self.columnCount = pheromoneDimension[1]
    return pheromone


def printPeptideMap(highestScoredPeptideMap):
    for k,v in highestScoredPeptideMap.items():
        print(str(k) + "-> [",end="")
        for j in v:
            print(j.peptide + "(" + str(j.mass) + "->" + str(j.score) + ")",end="")
            print(",")
        print("]")

def printEvaluationCounters(AcoLocalFitnessEvaluation, fe):
    print("AcoLocalFitnessEvaluation = " + str(AcoLocalFitnessEvaluation))
    print("fitnessEvaluationCount = " + str(fe))
