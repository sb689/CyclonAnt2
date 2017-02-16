from copy import copy,deepcopy
import numpy as np
import datetime
import sys
from Init import Init
from PeptideSequence import PeptideSequence
from BatchTest import BatchTest

AcoLocalFitnessEvaluation = 0
LogInterval = 8000

def exTestScript():
    global AcoLocalFitnessEvaluation
    global LogInterval
    batchTestObj = BatchTest()
    originalPeptides = batchTestObj.getOriginalPeptides()
    resultString = ""

    for i in range(3):
        testSpectrums = batchTestObj.getTestSpectrums(i)
        for j in range(len(testSpectrums)):
            spectrumIn = testSpectrums[j]
            print("length = " + str(len(spectrumIn)))
            realPeptide = originalPeptides[j]
            params = batchTestObj.selectParameters(len(spectrumIn))
            '''format of params -
            params[0] = acoIteration
            params[1] = finalLSIteration
            params[2] = M
            params[3] = N
            params[4] = ants
            params[5] = testCounter'''
            initObj = Init(spectrumIn,params[2],params[3],params[4])
            batchTestObj.setScoreObject(initObj)
            resultList = getResult(initObj,1,params[0],params[1],params[5])
            subString = batchTestObj.processResult(resultList, realPeptide)
            resultString = resultString + subString
            print("subString = " + subString)
            print(resultString)

            resultList = getResult(initObj,2,params[0],params[1],params[5])
            subString = batchTestObj.processResult(resultList, realPeptide)
            resultString = resultString + subString
            print("subString = " + subString)
            print(resultString)

            resultList = getResult(initObj,3,params[0],params[1],params[5])
            subString = batchTestObj.processResult(resultList, realPeptide)
            resultString = resultString + subString
            print("subString = " + subString)
            print(resultString)

            resultList = getResult(initObj,4,params[0],params[1],params[5])
            subString = batchTestObj.processResult(resultList, realPeptide)
            resultString = resultString + subString + "\n"
            print("subString = " + subString)
            print(resultString)

        resultString = resultString + "\n"
    print("resultString = ")
    print(resultString)



def getResult(initObj,method,acoIteration,finalLSIteration,testCounter):
    global AcoLocalFitnessEvaluation
    global LogInterval
    resultList = []
    for j in range(testCounter):
        peptideSequenceObj = PeptideSequence(initObj,method,acoIteration,finalLSIteration,
        AcoLocalFitnessEvaluation,LogInterval)
        result = deepcopy(peptideSequenceObj.executeMethod())
        resultList.append(result)
        fsock.flush()
    if method == 1:
        AcoLocalFitnessEvaluation = peptideSequenceObj.getAcoLocalFitnessEvaluation()

    return resultList


def main():
    global AcoLocalFitnessEvaluation
    global LogInterval
    """Main entry point for the script."""
    M,N,ants,Reads = sys.stdin.read().splitlines()
    M = int(M)
    N = int(N)
    ants = int(ants)
    Reads = Reads.split(",")
    spectrum = []
    spectrum = [int(x) for x in Reads]
    start = datetime.datetime.now()
    initObj = Init(spectrum,M,N,ants)
    peptideSequenceObj = PeptideSequence(initObj,1,1000,100,AcoLocalFitnessEvaluation,LogInterval)
    peptideSequenceObj.executeMethod()
    fsock.flush()
    AcoLocalFitnessEvaluation = peptideSequenceObj.getAcoLocalFitnessEvaluation()
    peptideSequenceObj = PeptideSequence(initObj,2,1000,100,AcoLocalFitnessEvaluation,LogInterval)
    peptideSequenceObj.executeMethod()
    fsock.flush()
    peptideSequenceObj = PeptideSequence(initObj,3,1000,100,AcoLocalFitnessEvaluation,LogInterval)
    peptideSequenceObj.executeMethod()
    fsock.flush()
    peptideSequenceObj = PeptideSequence(initObj,4,1000,100,AcoLocalFitnessEvaluation,LogInterval)
    peptideSequenceObj.executeMethod()
    fsock.flush()

    end = datetime.datetime.now()
    print("total time required = " + str(end-start))


if __name__ == '__main__':
    fsock = open("E:\\UIU_MSC\\bioInformatics\\GenomeAssembly\\aco_peptide\\output44.txt", 'w')
    saveout = sys.stdout
    sys.stdout = fsock
    main()
    #exTestScript()
    sys.stdout = saveout
    fsock.close()
    sys.exit()
#main segment ends
