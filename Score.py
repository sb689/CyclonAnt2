import numpy as np
from copy import copy,deepcopy

class Score(object):
    def __init__(self, spectrumMapIn):
        self.spectrumMap = spectrumMapIn


    def sequenceScoreExt(self,peptide):
        score = 0
        sequenceBreakDownMap = {}
        subsequenceMap = {}
        self.sequenceSpectrum(peptide, sequenceBreakDownMap, subsequenceMap)
        tmpSpectrumMap = copy(self.spectrumMap)
        keys = list(subsequenceMap.keys())
        for i in range(len(keys)):
            subsequence = keys[i]
            subsequenceList = subsequenceMap[keys[i]]
            if subsequence in tmpSpectrumMap:
                spectrumMapValues = tmpSpectrumMap[subsequence]
                if spectrumMapValues > 0:
                    score += 1
                else:
                    score = score-1 if score > 0 else score
                tmpSpectrumMap[subsequence] = spectrumMapValues - 1
            else:
                if score > 0 :
                    score = score-1 if score > 0 else score


            #no of mismatch
            mismatchNo = 0
            for j in range(len(subsequenceList)):

                if subsequenceList[j] in tmpSpectrumMap:
                    spectrumMapValues = tmpSpectrumMap[subsequenceList[j]]
                    if spectrumMapValues > 0:
                        score += 2
                    else:
                        score = score-2 if score > 1 else score
                        mismatchNo += 1
                    tmpSpectrumMap[subsequenceList[j]] = spectrumMapValues - 1
                else:
                    score = score-2 if score > 1 else score
                    mismatchNo += 1

            if mismatchNo > int(len(subsequenceList)*1/2):
                score = score - 5
        return score



    def sequenceSpectrum(self,peptide,sequenceBreakDownMap,subsequenceMap):
        #print("sequenceSpectrum testing")
        prefixMass = []
        peptides = []
        peptides = peptide.split("-")
        peptides = [int(x) for x in peptides]
        n = len(peptides)
        for i in range(n):
            tmp = peptides[i]
            tmpList = []
            tmpList.append(peptides[i])
            for j in range(i+1,n):
                tmp += peptides[j]
                tmpList.append(peptides[j])
                if tmp not in sequenceBreakDownMap:
                    sequenceBreakDownMap[tmp] = []
                sequenceBreakDownMap[tmp].extend(tmpList)
        for i in range(n):
            tmp = peptides[i]
            if tmp not in subsequenceMap:
                subsequenceMap[tmp] = []
            tmpList = []
            for j in range(i+1,n):
                tmp += peptides[j]
                tmpList.append(tmp)
            subsequenceMap[peptides[i]].extend(tmpList)
        if 0 not in subsequenceMap:
            subsequenceMap[0] = []



    def sequenceScore(self,peptide):
        score = 0
        sequenceBreakDownMap = {}
        subsequenceMap = {}
        self.sequenceSpectrum(peptide, sequenceBreakDownMap, subsequenceMap)
        keys = list(subsequenceMap.keys())
        for i in range(len(keys)):
            subsequence = keys[i]
            subsequenceList = subsequenceMap[keys[i]]
            if subsequence in self.spectrumMap:
                score += 1
                for j in range(len(subsequenceList)):
                    if subsequenceList[j] in self.spectrumMap:
                        score += 2

            #print("after subsequence = " + str(subsequence) + " score = " + str(score))
        #print("sequence score = " + str(score))
        return score


    def linearScore(self,peptide):
        score = 0
        testSpectrum = linearSpectrum(peptide)
        #print("linearSpectrum")
        #print(testSpectrum)
        prevSubpeptide = ""
        for subpeptide in testSpectrum:
            if subpeptide == prevSubpeptide:
                continue
            subpeptideCount = testSpectrum.count(subpeptide)
            try:
                tmp = self.spectrumMap[subpeptide]
            except Exception:
                tmp = 0
            if  tmp <= subpeptideCount and tmp > 0:
                score += tmp
            elif tmp > subpeptideCount:
                score += subpeptideCount
            prevSubpeptide = subpeptide
        return score


    def cycloScore(self,peptide):
        score = 0
        if peptide == "":
            return 0
        testSpectrum = cyclicSpectrum(peptide)
        prevSubpeptide = ""

        for subpeptide in testSpectrum:
            if subpeptide == prevSubpeptide:
                continue
            subpeptideCount = testSpectrum.count(subpeptide)
            #print("for subpeptide " + str(subpeptide) + "subpeptideCount is " + str(subpeptideCount))
            try:
                tmp = self.spectrumMap[subpeptide]
            except Exception:
                tmp = 0
            if  tmp <= subpeptideCount and tmp > 0:
                score += tmp
            elif tmp > subpeptideCount:
                score += subpeptideCount
            #print("tmp = " + str(tmp))
            #print("score = " + str(score))
            prevSubpeptide = subpeptide
        return score


def linearSpectrum(peptide):
    linearSpectrum = []
    prefixMass = []
    peptides = []
    if type(peptide) is np.int32:
        peptides.append(peptide)
        n = 1
    else:
        peptides = peptide.split("-")
        n = len(peptides)
    for i in range(n+1):
        prefixMass.append(0)
    for i in range(n):
        mass = int(peptides[i])
        prefixMass[i+1] = prefixMass[i] + mass
    for i in range(n):
        for j in range(i+1,n+1):
            tmp = prefixMass[j] - prefixMass[i]
            linearSpectrum.append(tmp)
    linearSpectrum.append(0)
    linearSpectrum.sort()
    return linearSpectrum


def cyclicSpectrum(peptide):
    cyclicSpectrum = []
    prefixMass = []
    peptides = []
    if type(peptide) is np.int32:
        peptides.append(peptide)
        n = 1
    else:
        peptides = peptide.split("-")
        n = len(peptides)
    for i in range(n+1):
        prefixMass.append(0)
    for i in range(n):
        mass = int(peptides[i])
        prefixMass[i+1] = prefixMass[i] + mass
    peptideMass = prefixMass[n]
    for i in range(n):
        for j in range(i+1,n+1):
            tmp = prefixMass[j] - prefixMass[i]
            cyclicSpectrum.append(tmp)
            if i > 0 and j < n :
                tmp2 = peptideMass - tmp
                cyclicSpectrum.append(tmp2)
    cyclicSpectrum.append(0)
    cyclicSpectrum.sort()
    return cyclicSpectrum
