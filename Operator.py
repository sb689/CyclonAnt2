from copy import copy,deepcopy
import numpy as np
from FitnessEvaluation import FitnessEvaluation
from Score import Score
from PeptideItem import PeptideItem
from Utility import Utility

class Operator(object):
    def __init__(self, fitnessObjIn, scoreObjIn):
        self.fitnessObj = fitnessObjIn
        self.scoreObj = scoreObjIn
        self.utilityObj = Utility()
        self.swapSuccess = 0
        self.fissonSuccess = 0
        self.fusionSuccess = 0
        self.transportationSuccess = 0
        self.mergeAndBreakSuccess = 0


    def fusion(self, peptide, extensiveScoring, mapToConcatenate):
        splittedPeptide = copy(peptide.peptide.split("-"))
        splittedPeptide = [int(x) for x in splittedPeptide]
        concatKeys = list(mapToConcatenate.keys())
        #print("concatKeys = ")
        #print(concatKeys)
        isChanged = False
        resultlist = []
        comparisonMap = {}

        for firstFragmentIndex in range(len(splittedPeptide)):
            if splittedPeptide[firstFragmentIndex] in concatKeys:
                #print(str(splittedPeptide[firstFragmentIndex]) + " first frag is present in splittedPeptide")
                concatValues = mapToConcatenate[splittedPeptide[firstFragmentIndex]]
                #print("concatValues = ")
                #print(concatValues)
                for secondFragment in range(len(concatValues)):
                    if concatValues[secondFragment] in splittedPeptide:
                        #print(str(concatValues[secondFragment]) + " second frag is present in splittedPeptide")
                        summationOfFragments = splittedPeptide[firstFragmentIndex] + concatValues[secondFragment]
                        secondFragmentIndex = splittedPeptide.index(concatValues[secondFragment])
                        newPeptide = copy(splittedPeptide)
                        newPeptide[firstFragmentIndex] = summationOfFragments
                        del newPeptide[secondFragmentIndex]
                        generatedPeptide = PeptideItem()
                        for k in range(len(newPeptide)):
                            generatedPeptide.peptide += str(newPeptide[k]) + "-"
                        generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
                        if extensiveScoring:
                            generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
                        else:
                            generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)
                        self.fitnessObj.increaseFE(generatedPeptide)

                        if generatedPeptide.score > peptide.score:
                            self.fusionSuccess += 1
                            #print("fussion succeed yee")
                        #print("generated intermidiate peptide in Fusion")
                        #print(generatedPeptide.peptide + " score = " + str(generatedPeptide.score))
                            if generatedPeptide.score not in comparisonMap:
                                comparisonMap[generatedPeptide.score] = []
                            comparisonMap[generatedPeptide.score].append(generatedPeptide)
        #print("from Fusion")
        bestScoredPeptideList = []
        if len(comparisonMap) > 0:
            bestScoredPeptideList = getBestScoredPeptidesFromMap(comparisonMap,peptide.score,False)
            for k in range(len(bestScoredPeptideList)):
                bestScoredPeptideList[k].mass = self.utilityObj.calculateMassOfPeptide(bestScoredPeptideList[k].peptide)
        return bestScoredPeptideList


    def fission(self,peptide, extensiveScoring,mapToSplit):
        #print("incoming peptide in Fission")
        #print(peptide.peptide + " score = " + str(peptide.score))
        splittedPeptide = copy(peptide.peptide.split("-"))
        splittedPeptide = [int(x) for x in splittedPeptide]
        splitKeys = list(mapToSplit.keys())
        #print("splitKeys = ")
        #print(splitKeys)
        isChanged = False
        resultlist = []
        comparisonMap = {}

        for i in range(len(splittedPeptide)):
            if splittedPeptide[i] in splitKeys:
                splitValues = mapToSplit[splittedPeptide[i]]
                #print("for "+ str(splittedPeptide[i]) +" len of splitValues = " + str(len(splitValues)))
                #print(splitValues)
                for j in range(len(splitValues)):
                    valueToBeReplacedWith = splitValues[j]
                    newPeptide = copy(splittedPeptide)
                    del newPeptide[i]
                    for k in range(len(valueToBeReplacedWith)):
                        newPeptide.insert(i+k,valueToBeReplacedWith[k])
                    generatedPeptide = PeptideItem()
                    for k in range(len(newPeptide)):
                        generatedPeptide.peptide += str(newPeptide[k]) + "-"
                    generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
                    if extensiveScoring:
                        generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
                    else:
                        generatedPeptide.score = self.scoreObj.SequenceScore(generatedPeptide.peptide)
                    self.fitnessObj.increaseFE(generatedPeptide)

                    if generatedPeptide.score > peptide.score:
                        self.fissonSuccess += 1
                        #print("Fission succeed yoo")
                    #print("generated intermidiate peptide in Fission")
                    #print(generatedPeptide.peptide + " score = " + str(generatedPeptide.score))
                        if generatedPeptide.score not in comparisonMap:
                            comparisonMap[generatedPeptide.score] = []
                        comparisonMap[generatedPeptide.score].append(generatedPeptide)
        #print("from Fission")
        bestScoredPeptideList = []
        if len(comparisonMap) > 0:
            bestScoredPeptideList = getBestScoredPeptidesFromMap(comparisonMap,peptide.score,False)
            for k in range(len(bestScoredPeptideList)):
                bestScoredPeptideList[k].mass = self.utilityObj.calculateMassOfPeptide(bestScoredPeptideList[k].peptide)
                bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)
        return bestScoredPeptideList


    def swap(self,peptide, extensiveScoring):
        #print("incoming peptide in Swap, " + peptide.peptide + "(" +str(peptide.score) + ")")
        comparisonMap = {}
        for i in range(5):
            splittedPeptide = copy(peptide.peptide.split("-"))
            rand1 = np.random.randint(0,len(splittedPeptide))
            rand2 = np.random.randint(0,len(splittedPeptide))
            while splittedPeptide[rand1] == splittedPeptide[rand2]:
                rand2 = np.random.randint(0,len(splittedPeptide))
            holdAminoAcid = splittedPeptide[rand2]
            splittedPeptide[rand2] = splittedPeptide[rand1]
            splittedPeptide[rand1] = holdAminoAcid

            '''rand1 = np.random.randint(0,len(splittedPeptide)-1)
            holdAminoAcid = splittedPeptide[rand1+1]
            splittedPeptide[rand1+1] = splittedPeptide[rand1]
            splittedPeptide[rand1] = holdAminoAcid'''
            generatedPeptide = PeptideItem()
            for k in range(len(splittedPeptide)):
                generatedPeptide.peptide += str(splittedPeptide[k]) + "-"
            generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
            if extensiveScoring:
                generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
            else:
                generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)

            self.fitnessObj.increaseFE(generatedPeptide)
            if generatedPeptide.score > peptide.score:
                self.swapSuccess += 1
                if generatedPeptide.score not in comparisonMap:
                    comparisonMap[generatedPeptide.score] = []
                comparisonMap[generatedPeptide.score].append(generatedPeptide)

        #print("from Swap")
        bestScoredPeptideList = []
        if len(comparisonMap) > 0:
            bestScoredPeptideList = getBestScoredPeptidesFromMap(comparisonMap,peptide.score,True)
            for k in range(len(bestScoredPeptideList)):
                bestScoredPeptideList[k].mass = self.utilityObj.calculateMassOfPeptide(bestScoredPeptideList[k].peptide)
                bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)
        return bestScoredPeptideList


    def transposition(self,peptide, extensiveScoring):
        #print("incoming peptide in Transposition, " + peptide.peptide + "(" +str(peptide.score) + ")")
        comparisonMap = {}

        for i in range(5):
            splittedPeptide = copy(peptide.peptide.split("-"))
            #print(splittedPeptide)
            rand2 = np.random.randint(0,len(splittedPeptide))
            holdAminoAcid = splittedPeptide[rand2]
            del splittedPeptide[rand2]
            rand1 = np.random.randint(0,len(splittedPeptide))
            while holdAminoAcid == splittedPeptide[rand1]:
                rand1 = np.random.randint(0,len(splittedPeptide))
            splittedPeptide.insert(rand1,holdAminoAcid)

            generatedPeptide = PeptideItem()
            for k in range(len(splittedPeptide)):
                generatedPeptide.peptide += str(splittedPeptide[k]) + "-"
            generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
            if extensiveScoring:
                generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
            else:
                generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)
            self.fitnessObj.increaseFE(generatedPeptide)
            if generatedPeptide.score > peptide.score:
                self.transportationSuccess += 1
                if generatedPeptide.score not in comparisonMap:
                    comparisonMap[generatedPeptide.score] = []
                comparisonMap[generatedPeptide.score].append(generatedPeptide)

        #print("from Transposition")
        bestScoredPeptideList = []
        if len(comparisonMap) > 0:
            bestScoredPeptideList = getBestScoredPeptidesFromMap(comparisonMap,peptide.score,True)
            for k in range(len(bestScoredPeptideList)):
                bestScoredPeptideList[k].mass = self.utilityObj.calculateMassOfPeptide(bestScoredPeptideList[k].peptide)
                bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)
        return bestScoredPeptideList


    def mutate(self,peptide,extensiveScoring,initObj):
        subsequenceMap = {}
        sequenceBreakDownMap = {}
        self.scoreObj.sequenceSpectrum(peptide.peptide, sequenceBreakDownMap,subsequenceMap)
        comparisonMap = {}
        comparisonMap[peptide.score] = []
        comparisonMap[peptide.score].append(peptide)
        #print("incoming peptide in Mutate")
        #print(peptide.peptide + "(" + str(peptide.score) + ")")
        #print("sequenceBreakDownMap")
        #print(sequenceBreakDownMap)
        #print("subsequenceMap")
        #print(subsequenceMap)
        #print("peptide = " + peptide.peptide)
        #print("===================")
        aminoAcidToChange = -1
        keys = list(sequenceBreakDownMap.keys())
        keys.sort(reverse = True)
        shouldStop = False
        for i in range(len(keys)):
            isKeyPresent = False
            #print("keys[i] = " + str(keys[i]))
            if keys[i] in initObj.spectrumMap:
                isKeyPresent = True
            subsequenceList = sequenceBreakDownMap[keys[i]]
            #print("subsequenceList")
            #print(subsequenceList)
            for j in range(len(subsequenceList)):
                if subsequenceList[j] not in initObj.spectrumMap:
                    sequenceList = subsequenceMap[subsequenceList[j]]
                    #print("sequenceList for " + str(subsequenceList[j]))
                    #print(sequenceList)
                    noOfSequencesPresent = 0
                    for k in range(len(sequenceList)):
                        if sequenceList[k] in initObj.spectrumMap:
                            noOfSequencesPresent += 1
                    #print("noOfSequencesPresent for " + str(subsequenceList[j])  +" is = " + str(noOfSequencesPresent))
                    if noOfSequencesPresent < int(len(sequenceList)*2/3):
                        aminoAcidToChange = subsequenceList[j]
                        shouldStop = True
                        break
            if shouldStop:
                break
        splittedPeptide = copy(peptide.peptide.split("-"))
        splittedPeptide = [int(x) for x in splittedPeptide]
        #print("aminoAcidToChange = " + str(aminoAcidToChange))
        if aminoAcidToChange != -1:
            positionTochange = splittedPeptide.index(aminoAcidToChange)
        else:
            positionTochange = np.random.randint(0,len(splittedPeptide))

        for j in range(3):
            rand1 = np.random.randint(0,len(initObj.aminoAcidMass))
            splittedPeptide[positionTochange] = initObj.aminoAcidMass[rand1]
            generatedPeptide = PeptideItem()
            for k in range(len(splittedPeptide)):
                generatedPeptide.peptide += str(splittedPeptide[k]) + "-"
            generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
            if extensiveScoring:
                generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
            else:
                generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)

            #("generatedPeptide.score = " + str(generatedPeptide.score) + "," + "generatedPeptide = " + generatedPeptide.peptide)
            self.fitnessObj.increaseFE(generatedPeptide)
            if generatedPeptide.score not in comparisonMap:
                comparisonMap[generatedPeptide.score] = []
            comparisonMap[generatedPeptide.score].append(generatedPeptide)
        #print("from mutate")
        bestScoredPeptideList = getBestScoredPeptidesFromMap(comparisonMap,peptide.score,True)
        for k in range(len(bestScoredPeptideList)):
            bestScoredPeptideList[k].mass = self.utilityObj.calculateMassOfPeptide(bestScoredPeptideList[k].peptide)
            bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)
        return bestScoredPeptideList


    def reverse(self,peptideItemList,extensiveScoring):
        reversePeptideList = []
        for i in range(len(peptideItemList)):
            peptide = peptideItemList[i]
            splittedPeptide = peptide.peptide.split("-")
            splittedPeptide = [int(x) for x in splittedPeptide]
            reversedpeptide = ""
            for k in range(len(splittedPeptide)):
                reversedpeptide = str(splittedPeptide[k]) + "-" + reversedpeptide
            reversedpeptide = reversedpeptide.strip("-")
            generatedPeptide = PeptideItem()
            generatedPeptide.peptide = reversedpeptide
            if extensiveScoring:
                generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
            else:
                generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)
            generatedPeptide.mass = self.utilityObj.calculateMassOfPeptide(generatedPeptide.peptide)
            reversePeptideList.append(generatedPeptide)

        return reversePeptideList



    def mergeAndBreak(self,peptide,extensiveScoring,mergeAndBreakMap):
        #print("incoming peptide in MergeAndBreak, " + peptide.peptide + "(" +str(peptide.score) + ")")
        comparisonMap = {}
        splittedPeptide = copy(peptide.peptide.split("-"))
        splittedPeptide = [int(x) for x in splittedPeptide]
        mapKeys = list(mergeAndBreakMap.keys())
        shouldBreak = False
        for i in range(len(splittedPeptide)-1):
            sumOfMass = splittedPeptide[i] + splittedPeptide[i+1]
            #print("sumOfMass = " + str(sumOfMass))
            if sumOfMass in mapKeys:
                #print("sumOfMass found in mapKeys")
                mapValues = mergeAndBreakMap[sumOfMass]
                for j in range(len(mapValues)):
                    pair = mapValues[j]
                    if splittedPeptide[i] not in pair:
                        splittedPeptide2 = copy(splittedPeptide)
                        splittedPeptide2[i] = pair[0]
                        splittedPeptide2[i+1] = pair[1]
                        generatedPeptide = PeptideItem()
                        for k in range(len(splittedPeptide2)):
                            generatedPeptide.peptide += str(splittedPeptide2[k]) + "-"
                        generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
                        if extensiveScoring:
                            generatedPeptide.score = self.scoreObj.sequenceScoreExt(generatedPeptide.peptide)
                        else:
                            generatedPeptide.score = self.scoreObj.sequenceScore(generatedPeptide.peptide)
                        #print("generatedPeptide in MergeAndBreak")
                        #print(generatedPeptide.peptide + "->(" + str(generatedPeptide.score) + ")")
                        self.fitnessObj.increaseFE(generatedPeptide)
                        if generatedPeptide.score > peptide.score:
                            generatedPeptide.mass = self.utilityObj.calculateMassOfPeptide(generatedPeptide.peptide)
                            shouldBreak = True
                            self.mergeAndBreakSuccess += 1
                            break
            if shouldBreak:
                break
        bestScoredPeptideList = []
        if shouldBreak:
            bestScoredPeptideList.append(generatedPeptide)
        return bestScoredPeptideList


def getBestScoredPeptidesFromMap(comparisonMap,peptideScore,considerEqual):
    keys = list(comparisonMap.keys())
    keys.sort(reverse = True)
    '''print("comparisonMap")
    for k in keys:
        resultList = comparisonMap[k]
        print(str(k) + "->[",end="")
        for r in resultList:
            print(r.peptide + ",")
        print("]")'''
    resultlist = []
    if len(keys)>0:
        maxKey = keys[0]
        if considerEqual:
            if maxKey >= peptideScore:
                resultlist = comparisonMap[maxKey]
                resultlist = list(set(resultlist))
        else:
            if maxKey > peptideScore:
                resultlist = comparisonMap[maxKey]
                resultlist = list(set(resultlist))

    return resultlist
