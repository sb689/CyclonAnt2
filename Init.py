import numpy as np

class Init(object):
    def __init__(self,spectrumIn,M,N,antsIn):
        self.m = M
        self.n = N
        self.ants = antsIn
        self.spectrum = spectrumIn
        self.mergeAndBreakMap = {}
        self.mapToSplit = {}
        self.mapToConcatenate = {}
        self.spectrumMap = {}
        self.aminoAcidMass = []
        self.parentMass = 0
        self.pheromoneConstant = 0.7
        self.beta = 2
        self.q = 0.1


        self.parentMass = max(self.spectrum)
        self.spectrumMap = populateSpectrumMap(self.spectrum)
        self.aminoAcidMass = getFrequentAminoAcids(self.spectrum, self.m)
        self.mergeAndBreakMap = createAminoAcidMergeAndBreakMap(self.aminoAcidMass)
        createAminoAcidRearrangementMap(self.aminoAcidMass, self.mapToConcatenate, self.mapToSplit)
        '''print("spectrumMap")
        print(self.spectrumMap)
        print("aminoAcidMass")
        print(self.aminoAcidMass)
        print("mergeAndBreakMap")
        print(self.mergeAndBreakMap)
        print("mapToConcatenate")
        print(self.mapToConcatenate)
        print("mapToSplit")
        print(self.mapToSplit)'''


'''def populateSpectrum(spectrumIn):
    spectrum = []
    spectrum = [int(x) for x in spectrumIn]
    return spectrum'''


def populateSpectrumMap(spectrum):
    spectrumMap = {}
    for i in range(len(spectrum)):
        item = spectrum[i]
        if item not in spectrumMap:
            spectrumMap[item] = 0
        spectrumMap[item] += 1
    #print("SpectrumMap")
    #print(SpectrumMap)
    return spectrumMap


def getFrequentAminoAcids(spectrum,m):
    output = []
    aminoAcidMass = []
    for i in range(len(spectrum)):
        for j in range(i):
            diff = abs(spectrum[j] - spectrum[i])
            if diff >= 57 and diff <=200:
                output.append(diff)
    #print("Length of Spectrum " + str(len(Spectrum)))
    #print("length of output " + str(outputLen))

    output.sort()
    prevMass = 0
    countMap = {}
    outputUnique = np.unique(output)
    for i in range(len(outputUnique)):
        mass = outputUnique[i]
        massCount = output.count(mass)
        if massCount not in countMap:
            countMap[massCount] = []
        countMap[massCount].append(mass)
    #print("countMap , len "+ str(len(countMap)))
    #print('\n'.join(str(k)+ " -> " + str(v) for k,v in countMap.items()))

    keys = list(countMap.keys())
    keys.sort(reverse = True)
    for index in range(len(keys)):
        tmpArr = countMap[keys[index]]
        aminoAcidMass.extend(tmpArr)
        if len(aminoAcidMass) >= m:
            break
    aminoAcidMass.sort()
    return aminoAcidMass


def createAminoAcidMergeAndBreakMap(aminoAcidMass):
    mergeAndBreakMap = {}
    for i in range(len(aminoAcidMass)-1):
        firstAminoAcid = aminoAcidMass[i]
        for j in range(i+1,len(aminoAcidMass)):
            secondAminoAcid = aminoAcidMass[j]
            sumOfMass = firstAminoAcid + secondAminoAcid
            if sumOfMass not in mergeAndBreakMap:
                mergeAndBreakMap[sumOfMass] = []
            mergeAndBreakMap[sumOfMass].append([firstAminoAcid,secondAminoAcid])
    return mergeAndBreakMap


def createAminoAcidRearrangementMap(aminoAcidMass, mapToConcatenate, mapToSplit):
    for i in range(len(aminoAcidMass)):
        k = aminoAcidMass[i]
        for j in range(i+1,len(aminoAcidMass)):
            v = aminoAcidMass[j]
            sumOfMass = k + v
            if sumOfMass > 200:
                break
            if sumOfMass in aminoAcidMass:
                if k not in mapToConcatenate:
                    mapToConcatenate[k] = []
                mapToConcatenate[k].append(v)
                if v not in mapToConcatenate:
                    mapToConcatenate[v] = []
                mapToConcatenate[v].append(k)
                if sumOfMass not in mapToSplit:
                    mapToSplit[sumOfMass] = []
                mapToSplit[sumOfMass].append([k,v])
