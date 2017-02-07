from copy import copy,deepcopy
import numpy as np
import datetime
import sys

#this version compares algorithm with updated fitness evaluation. but the compared score is CycloScore.
#also it takes average score from 5 runs(not in Batchtest, the other way around)

MergeAndBreakMap={}
MapToSplit = {}
MapToConcatenate = {}
AminoAcidMass = []
SpectrumMap = {}
PheromoneConstant = 0.7
Beta = 2
Q = 0.1
LocalSearchIteration = 1
Pheromone = []
ParentMass = 0
HighestScoredPeptideMap = {}
Spectrum = []
SwapSuccess = 0
FissonSuccess = 0
FussionSuccess = 0
TransportationSuccess = 0
MergeAndBreakSuccess = 0
AcoLocalFitnessEvaluation = 0
LogInterval = 8000
score1 = []
score2 = []
score3 = []


class PeptideItem:
    def __init__(self):
        self.mass = 0
        self.peptide = ""
        self.score = 0


    def printPeptide(self):
        print("peptide = " + self.peptide)
        print("mass = " + str(self.mass))
        print("score = " + str(self.score))


class FitnessEvaluation():
    def __init__(self,method):
        global LogInterval
        global score1
        global score2
        global score3
        self.fitnessEvaluationCount = 0
        self.scoreStr = ""
        self.fitnessStr = ""
        self.method = method
        self.nextLogCounter = LogInterval
        self.maxScoredPeptide = PeptideItem()
        self.partialScore = []

    def increaseFE(self, newPeptide):
        self.fitnessEvaluationCount += 1
        if newPeptide.score > self.maxScoredPeptide.score:
            self.maxScoredPeptide = copy(newPeptide)

        if self.fitnessEvaluationCount >= self.nextLogCounter :
            #cyScore = CycloScore(self.maxScoredPeptide.peptide)
            self.scoreStr = self.scoreStr + str(self.maxScoredPeptide.score) + ";"
            self.fitnessStr = self.fitnessStr + str(self.fitnessEvaluationCount) + ";"
            #self.partialScore.append(self.maxScoredPeptide.score)
            self.nextLogCounter += LogInterval



    def getFitnessEvaluationCount(self):
        return self.fitnessEvaluationCount


    def printFE(self):
        print("FitnessEvaluation count  for method = " + str(self.method))
        print(self.scoreStr)
        print(self.fitnessStr)
        '''if self.method == 1:
            score1.append(self.partialScore)
        elif self.method == 2:
            score2.append(self.partialScore)
        else:
            score3.append(self.partialScore)'''



#Batch test segment starts here

def GetOriginalPeptides():
    originalPeptides = [0,0,0,0,0,0]
    originalPeptides[0] = "97-163-99-97-113-113-113-97"
    originalPeptides[1] = "99-114-113-147-97-99-114-113-147-97"
    originalPeptides[2] = "99-114-113-147-97-147-147-114-128-163"
    originalPeptides[3] = "99-114-113-147-97-163-99-114-113-147-97-163"
    originalPeptides[4] = "97-115-129-163-115-163-87-115-113-129-115-71-115"
    originalPeptides[5] = "99-114-113-147-97-163-99-114-113-147-97-163-88-96-174"
    return originalPeptides


def GetTestSpectrums(tag):
    #1,2,3
    originalSpectrum = [0,0,0,0,0,0]
    spectrum_10 = [0,0,0,0,0,0]
    spectrum_25 = [0,0,0,0,0,0]

    if tag == 0:
        originalSpectrum[0] = [0, 97, 97, 97, 99, 113, 113, 113, 163, 194, 196, 210, 210, 226, 226, 260, 262, 307, 309, 323, 323, 339, 357, 359, 359, 420, 422, 436, 436, 456, 456, 470, 472, 533, 533, 535, 553, 569, 569, 583, 585, 630, 632, 666, 666, 682, 682, 696, 698, 729, 779, 779, 779, 793, 795, 795, 795, 892]
        originalSpectrum[1] = [0, 97, 97, 99, 99, 113, 113, 114, 114, 147, 147, 196, 196, 213, 213, 227, 227, 244, 244, 260, 260, 310, 310, 326, 326, 343, 343, 357, 357, 374, 374, 423, 423, 456, 456, 457, 457, 471, 471, 473, 473, 570, 570, 570, 570, 570, 570, 570, 570, 570, 570, 667, 667, 669, 669, 683, 683, 684, 684, 717, 717, 766, 766, 783, 783, 797, 797, 814, 814, 830, 830, 880, 880, 896, 896, 913, 913, 927, 927, 944, 944, 993, 993, 1026, 1026, 1027, 1027, 1041, 1041, 1043, 1043, 1140]
        originalSpectrum[2] = [0, 97, 99, 113, 114, 114, 128, 147, 147, 147, 163, 213, 227, 242, 244, 244, 260, 261, 262, 291, 294, 326, 357, 374, 376, 389, 390, 391, 391, 405, 408, 471, 473, 489, 504, 504, 504, 505, 536, 538, 552, 570, 617, 618, 618, 633, 636, 651, 651, 652, 699, 717, 731, 733, 764, 765, 765, 765, 780, 796, 798, 861, 864, 878, 878, 879, 880, 893, 895, 912, 943, 975, 978, 1007, 1008, 1009, 1025, 1025, 1027, 1042, 1056, 1106, 1122, 1122, 1122, 1141, 1155, 1155, 1156, 1170, 1172, 1269]
        originalSpectrum[3] = [0, 97, 97, 99, 99, 113, 113, 114, 114, 147, 147, 163, 163, 213, 213, 227, 227, 244, 244, 260, 260, 260, 260, 262, 262, 326, 326, 357, 357, 359, 359, 374, 374, 376, 376, 407, 407, 471, 471, 473, 473, 473, 473, 489, 489, 506, 506, 520, 520, 570, 570, 586, 586, 619, 619, 620, 620, 634, 634, 636, 636, 733, 733, 733, 733, 733, 733, 733, 733, 733, 733, 733, 733, 830, 830, 832, 832, 846, 846, 847, 847, 880, 880, 896, 896, 946, 946, 960, 960, 977, 977, 993, 993, 993, 993, 995, 995, 1059, 1059, 1090, 1090, 1092, 1092, 1107, 1107, 1109, 1109, 1140, 1140, 1204, 1204, 1206, 1206, 1206, 1206, 1222, 1222, 1239, 1239, 1253, 1253, 1303, 1303, 1319, 1319, 1352, 1352, 1353, 1353, 1367, 1367, 1369, 1369, 1466]
        originalSpectrum[4] = [0, 71, 87, 97, 113, 115, 115, 115, 115, 115, 129, 129, 163, 163, 186, 186, 202, 212, 212, 228, 242, 244, 244, 250, 278, 278, 283, 292, 301, 315, 315, 327, 341, 357, 357, 365, 365, 398, 398, 407, 407, 428, 430, 441, 444, 456, 472, 478, 480, 504, 513, 522, 527, 527, 528, 543, 543, 559, 570, 593, 607, 619, 619, 630, 640, 642, 642, 643, 657, 658, 685, 690, 722, 722, 734, 745, 755, 755, 756, 771, 772, 772, 782, 793, 805, 805, 837, 842, 869, 870, 884, 885, 885, 887, 897, 908, 908, 920, 934, 957, 968, 984, 984, 999, 1000, 1000, 1005, 1014, 1023, 1047, 1049, 1055, 1071, 1083, 1086, 1097, 1099, 1120, 1120, 1129, 1129, 1162, 1162, 1170, 1170, 1186, 1200, 1212, 1212, 1226, 1235, 1244, 1249, 1249, 1277, 1283, 1283, 1285, 1299, 1315, 1315, 1325, 1341, 1341, 1364, 1364, 1398, 1398, 1412, 1412, 1412, 1412, 1412, 1414, 1430, 1440, 1456, 1527]
        originalSpectrum[5] = [0, 88, 96, 97, 97, 99, 99, 113, 113, 114, 114, 147, 147, 163, 163, 174, 184, 213, 213, 227, 227, 244, 244, 251, 260, 260, 260, 260, 262, 270, 273, 326, 326, 347, 348, 357, 357, 358, 359, 369, 374, 374, 376, 387, 407, 407, 444, 457, 471, 471, 473, 473, 473, 483, 489, 495, 500, 506, 520, 520, 521, 570, 570, 571, 586, 591, 596, 608, 618, 619, 620, 620, 634, 634, 636, 647, 684, 704, 717, 722, 733, 733, 733, 733, 733, 733, 733, 734, 743, 744, 765, 818, 821, 830, 831, 831, 832, 840, 846, 847, 847, 864, 878, 880, 896, 907, 917, 928, 944, 946, 960, 977, 977, 978, 984, 992, 993, 993, 994, 1003, 1006, 1059, 1080, 1081, 1090, 1091, 1091, 1091, 1091, 1091, 1091, 1091, 1102, 1107, 1120, 1140, 1177, 1188, 1190, 1190, 1204, 1204, 1205, 1206, 1216, 1228, 1233, 1238, 1253, 1254, 1254, 1303, 1304, 1304, 1318, 1324, 1329, 1335, 1341, 1351, 1351, 1351, 1353, 1353, 1367, 1380, 1417, 1417, 1437, 1448, 1450, 1450, 1455, 1465, 1466, 1467, 1467, 1476, 1477, 1498, 1498, 1551, 1554, 1562, 1564, 1564, 1564, 1564, 1573, 1580, 1580, 1597, 1597, 1611, 1611, 1640, 1650, 1661, 1661, 1677, 1677, 1710, 1710, 1711, 1711, 1725, 1725, 1727, 1727, 1728, 1736, 1824]
        print("returning originalSpectrum")
        return originalSpectrum
    elif tag == 1:
        spectrum_10[0] = [0, 97, 97, 97, 99, 113, 113, 113, 163, 194, 196, 210, 226, 226, 262, 307, 309, 323, 323, 357, 359, 359, 420, 422, 436, 436, 456, 456, 470, 472, 533, 533, 535, 569, 569, 583, 585, 630, 632, 666, 682, 682, 696, 698, 729, 779, 779, 779, 793, 795, 795, 795, 892]
        spectrum_10[1] = [0, 99, 99, 113, 113, 114, 114, 147, 147, 196, 196, 213, 213, 227, 244, 244, 260, 260, 310, 326, 326, 343, 357, 357, 374, 374, 423, 423, 456, 456, 457, 471, 471, 473, 473, 570, 570, 570, 570, 570, 570, 570, 570, 570, 667, 669, 669, 683, 683, 684, 684, 717, 717, 766, 766, 783, 783, 797, 797, 814, 814, 830, 830, 880, 896, 896, 913, 913, 927, 927, 944, 944, 993, 993, 1026, 1026, 1027, 1027, 1041, 1041, 1043, 1043, 1140]
        spectrum_10[2] = [0, 97, 99, 113, 114, 114, 128, 147, 163, 213, 227, 242, 244, 244, 260, 261, 262, 291, 326, 357, 374, 376, 389, 390, 391, 391, 405, 408, 471, 489, 504, 504, 504, 505, 536, 538, 552, 617, 618, 618, 633, 651, 651, 652, 699, 717, 731, 733, 764, 765, 765, 765, 796, 798, 861, 878, 878, 879, 880, 893, 895, 912, 943, 975, 978, 1007, 1008, 1009, 1025, 1025, 1027, 1042, 1056, 1106, 1122, 1122, 1141, 1155, 1155, 1156, 1170, 1172, 1269]
        spectrum_10[3] = [0, 97, 97, 99, 99, 113, 113, 114, 147, 147, 163, 213, 213, 227, 227, 244, 260, 260, 260, 260, 262, 326, 326, 357, 357, 359, 359, 374, 374, 376, 376, 407, 407, 471, 471, 473, 473, 473, 473, 489, 489, 506, 520, 570, 570, 586, 619, 619, 620, 620, 634, 634, 636, 636, 733, 733, 733, 733, 733, 733, 733, 733, 733, 733, 733, 830, 832, 832, 846, 846, 847, 847, 880, 880, 896, 896, 946, 946, 960, 960, 977, 977, 993, 993, 993, 993, 995, 1059, 1059, 1090, 1090, 1092, 1092, 1107, 1109, 1109, 1140, 1140, 1204, 1204, 1206, 1206, 1206, 1206, 1222, 1222, 1239, 1253, 1253, 1303, 1303, 1319, 1352, 1352, 1353, 1353, 1367, 1367, 1369, 1369, 1466]
        spectrum_10[4] = [0, 71, 87, 97, 113, 115, 115, 115, 115, 129, 129, 163, 163, 186, 186, 202, 212, 212, 228, 242, 244, 244, 250, 278, 278, 283, 292, 301, 315, 315, 327, 341, 357, 357, 365, 365, 398, 407, 407, 428, 430, 441, 444, 456, 472, 478, 480, 504, 513, 522, 527, 527, 528, 543, 559, 570, 593, 607, 619, 619, 630, 640, 642, 642, 643, 657, 685, 690, 722, 722, 734, 755, 755, 771, 772, 772, 782, 805, 805, 837, 842, 870, 884, 885, 885, 887, 897, 908, 908, 920, 934, 957, 968, 984, 984, 999, 1000, 1000, 1005, 1023, 1047, 1049, 1055, 1071, 1083, 1086, 1097, 1099, 1120, 1129, 1129, 1162, 1162, 1170, 1170, 1186, 1200, 1212, 1212, 1235, 1249, 1249, 1277, 1283, 1283, 1285, 1315, 1315, 1325, 1341, 1341, 1364, 1398, 1398, 1412, 1412, 1412, 1412, 1414, 1430, 1440, 1456, 1527]
        spectrum_10[5] = [0, 88, 96, 97, 97, 99, 99, 113, 113, 114, 114, 147, 147, 163, 184, 213, 227, 227, 244, 244, 251, 260, 260, 260, 260, 270, 326, 326, 347, 348, 357, 357, 358, 359, 369, 374, 374, 376, 387, 407, 457, 471, 471, 473, 473, 489, 495, 500, 506, 520, 520, 521, 570, 570, 571, 586, 591, 596, 608, 618, 619, 620, 620, 634, 634, 636, 647, 704, 717, 722, 733, 733, 733, 733, 733, 733, 733, 734, 743, 744, 765, 818, 821, 830, 831, 831, 832, 840, 846, 847, 847, 864, 878, 880, 896, 907, 917, 928, 944, 946, 977, 977, 984, 992, 993, 993, 994, 1003, 1006, 1059, 1080, 1081, 1090, 1091, 1091, 1091, 1091, 1091, 1091, 1091, 1102, 1107, 1120, 1140, 1177, 1188, 1190, 1204, 1204, 1205, 1206, 1216, 1228, 1233, 1238, 1253, 1254, 1254, 1303, 1304, 1304, 1318, 1324, 1329, 1335, 1341, 1351, 1351, 1351, 1353, 1353, 1367, 1380, 1417, 1437, 1448, 1450, 1455, 1465, 1466, 1467, 1476, 1477, 1498, 1551, 1554, 1562, 1564, 1564, 1564, 1564, 1573, 1580, 1611, 1611, 1640, 1650, 1661, 1661, 1677, 1677, 1710, 1711, 1711, 1725, 1725, 1727, 1727, 1728, 1736, 1824]
        print("returning spectrum_10")
        return spectrum_10
    else:
        spectrum_25[0] = [0, 97, 97, 97, 113, 113, 194, 196, 210, 210, 226, 226, 260, 262, 307, 309, 323, 339, 357, 359, 359, 422, 436, 436, 456, 470, 472, 533, 553, 569, 585, 632, 666, 666, 682, 696, 698, 729, 779, 779, 779, 795, 795, 892]
        spectrum_25[1] = [0, 99, 99, 113, 114, 114, 147, 147, 196, 196, 213, 213, 260, 310, 326, 326, 343, 357, 374, 374, 423, 423, 456, 456, 457, 471, 473, 473, 570, 570, 570, 570, 570, 570, 570, 570, 667, 669, 669, 683, 683, 684, 717, 766, 766, 783, 797, 797, 814, 830, 830, 880, 896, 896, 913, 927, 927, 944, 944, 993, 1026, 1026, 1027, 1027, 1041, 1041, 1043, 1043, 1140]
        spectrum_25[2] = [0, 97, 99, 114, 128, 147, 147, 213, 227, 242, 244, 244, 260, 261, 262, 291, 294, 326, 374, 376, 389, 391, 391, 405, 408, 471, 489, 504, 504, 504, 505, 536, 538, 552, 570, 617, 618, 618, 633, 636, 651, 717, 733, 764, 765, 780, 798, 861, 864, 878, 879, 880, 895, 912, 943, 975, 1007, 1008, 1009, 1027, 1056, 1106, 1122, 1122, 1155, 1156, 1170, 1172, 1269]
        spectrum_25[3] = [0, 99, 113, 114, 114, 147, 163, 213, 227, 227, 244, 244, 260, 260, 260, 262, 326, 357, 357, 359, 359, 374, 376, 407, 407, 471, 471, 473, 473, 489, 506, 506, 520, 570, 570, 586, 586, 620, 620, 634, 634, 636, 733, 733, 733, 733, 733, 733, 733, 733, 733, 830, 830, 832, 832, 846, 846, 847, 880, 880, 896, 896, 946, 960, 960, 977, 977, 993, 993, 993, 993, 995, 1059, 1090, 1090, 1092, 1107, 1107, 1109, 1109, 1140, 1140, 1204, 1204, 1206, 1206, 1206, 1222, 1222, 1253, 1303, 1319, 1319, 1352, 1352, 1353, 1367, 1367, 1369, 1369, 1466]
        spectrum_25[4] = [0, 71, 87, 97, 113, 115, 115, 115, 115, 129, 163, 163, 186, 202, 212, 212, 228, 242, 244, 244, 250, 278, 283, 301, 315, 327, 341, 357, 357, 365, 398, 398, 407, 407, 428, 430, 441, 444, 456, 472, 478, 513, 522, 527, 527, 528, 543, 543, 559, 570, 593, 607, 619, 619, 630, 640, 642, 642, 643, 657, 690, 722, 722, 734, 755, 755, 756, 771, 772, 772, 782, 805, 837, 842, 869, 884, 885, 885, 887, 897, 908, 908, 934, 968, 984, 984, 1000, 1005, 1014, 1047, 1055, 1071, 1083, 1086, 1099, 1120, 1129, 1129, 1162, 1186, 1212, 1235, 1244, 1249, 1283, 1285, 1315, 1325, 1341, 1364, 1364, 1398, 1412, 1412, 1414, 1430, 1440, 1456, 1527]
        spectrum_25[5] = [0, 88, 96, 97, 97, 99, 99, 113, 113, 114, 114, 147, 147, 163, 184, 213, 213, 227, 227, 244, 251, 260, 260, 326, 326, 347, 348, 357, 357, 358, 374, 374, 376, 387, 407, 407, 444, 457, 471, 473, 473, 483, 489, 495, 506, 520, 570, 586, 591, 596, 618, 619, 620, 620, 634, 634, 636, 647, 684, 704, 717, 722, 733, 733, 733, 733, 743, 744, 821, 830, 831, 831, 832, 846, 847, 847, 864, 878, 880, 896, 907, 928, 946, 960, 977, 977, 978, 984, 992, 993, 994, 1003, 1059, 1080, 1081, 1091, 1091, 1091, 1091, 1091, 1091, 1120, 1177, 1190, 1204, 1204, 1205, 1206, 1216, 1228, 1233, 1253, 1254, 1254, 1303, 1304, 1304, 1318, 1324, 1329, 1335, 1351, 1351, 1353, 1353, 1380, 1437, 1448, 1465, 1467, 1467, 1476, 1477, 1498, 1498, 1551, 1554, 1564, 1564, 1564, 1573, 1580, 1580, 1597, 1611, 1611, 1640, 1650, 1661, 1677, 1710, 1710, 1711, 1725, 1725, 1727, 1727, 1728, 1736, 1824]
        print("returning spectrum_25")
        return spectrum_25



def MatchPeptide(originalPeptide, generatedPeptide):
    isMatch = False
    splittedPep1 = originalPeptide.split("-")
    splittedPep1 = [int(x) for x in splittedPep1]
    splittedPep2 = generatedPeptide.split("-")
    splittedPep2 = [int(x) for x in splittedPep2]

    if len(splittedPep1) == len(splittedPep2):
        indexs = []
        for i in range(len(splittedPep2)):
            if splittedPep2[i] == splittedPep1[0]:
                indexs.append(i)
        for i in range(len(indexs)):
            tmpPep2 = []
            tmpPep2.extend(splittedPep2[indexs[i]:len(splittedPep2)])
            tmpPep2.extend(splittedPep2[0:indexs[i]])
            matchLength = 0
            for j in range(len(splittedPep1)):
                if splittedPep1[j] == tmpPep2[j]:
                    matchLength += 1
                else:
                    break
            if matchLength == len(splittedPep1):
                isMatch = True

    return isMatch


def ProcessResult(highestScoredMapList, realPeptide):
    global SpectrumMap
    returnString = ""
    scoreList = []
    for i in range(len(highestScoredMapList)):
        highestScoredMap = highestScoredMapList[i]
        if len(highestScoredMap) == 0:
            matched = 0
            returnString = returnString + str(-1) + ";" + str(-1) + ";" + str(-1) + ";" + str(-1) + ";" + str(matched) + ";"
            scoreList.append(0)
            continue
        keys = list(highestScoredMap.keys())
        keys.sort(reverse = True)
        keyCounter = 0
        limit = min(len(keys),5)
        while keyCounter < limit:
            maxScore = keys[keyCounter]
            maxScoredValues = highestScoredMap[maxScore]
            isMatch = 0
            for j in range(len(maxScoredValues)):
                generatedPeptide = maxScoredValues[j]
                isMatch = MatchPeptide(realPeptide, generatedPeptide.peptide)
                if isMatch:
                    break
            keyCounter += 1
            if isMatch:
                break
        peptideToBeScored = ""
        if isMatch:
            peptideToBeScored = maxScoredValues[j].peptide
        else:
            peptideToBeScored = maxScoredValues[0].peptide
        linearScore = LinearScore(peptideToBeScored)
        cycloScore = CycloScore(peptideToBeScored)
        sequenceScore = SequenceScore(peptideToBeScored)
        sequenceScoreExt = SequenceScoreExt(peptideToBeScored)
        matched = 0
        if isMatch:
            matched = 1
        scoreList.append(cycloScore)
        returnString = returnString + str(linearScore) + ";" + str(cycloScore) + ";" + str(sequenceScore) + ";" + str(sequenceScoreExt) + ";" + str(matched) + ";"
    avgScore = sum(scoreList)/len(scoreList)
    maxScore = max(scoreList)
    returnString = str(maxScore) + ";" + str(avgScore) + ";" + returnString
    return returnString



def SelectParameters(position):
    acoIteration = 1000
    finalLSIteration = 100
    M = 20
    N = 10
    if position == 0:
        acoIteration = 1000
        finalLSIteration = 100
    elif position in range(1,3):
        acoIteration = 1400
        finalLSIteration = 100
    elif position in range(3,5):
        acoIteration = 1600
        finalLSIteration = 200
    else:
        acoIteration = 2500
        finalLSIteration = 200
        M = 30
        N = 15
    result = []
    result.append(acoIteration)
    result.append(finalLSIteration)
    result.append(M)
    result.append(N)
    return result


def Batchtest():
    resultString = ""
    originalPeptides = GetOriginalPeptides()
    ants = 10
    testCounter = 5

    for specCounter in range(3):
        testSpectrumList = GetTestSpectrums(specCounter)
        for i in range(len(testSpectrumList)):
            global Spectrum
            scoreList = []
            averageSocre = 0
            maxScore = 0
            parameterList = SelectParameters(i)
            #first elemnt in parameterList is acoIteration values
            #second elemnt in parameterList is finalLSIteration
            #third elemnt in parameterList is M - number of masses in amino acid list
            #fourth elemnt in parameterList is N - number of peptides in final map
            acoIteration = parameterList[0]
            finalLSIteration = parameterList[1]
            M = parameterList[2]
            N = parameterList[3]
            realPeptide = originalPeptides[i]

            Spectrum = testSpectrumList[i]
            PopulateSpectrumMap()
            GetFrequentAminoAcids(M)
            CreateAminoAcidRearrangementMap()
            CreateAminoAcidMergeAndBreakMap()

            Spectrum = testSpectrumList[i]
            PopulateSpectrumMap()
            highestScoredMapList = ExSequencing(M,N,ants,acoIteration,finalLSIteration, 1, testCounter)
            subString = ProcessResult(highestScoredMapList, realPeptide)
            resultString = resultString + subString
            print("subString = " + subString)
            print(resultString)
            highestScoredMapList = ExSequencing(M,N,ants,acoIteration,finalLSIteration, 2, testCounter)
            subString = ProcessResult(highestScoredMapList, realPeptide)
            resultString = resultString + subString
            print("subString = " + subString)
            print(resultString)
            highestScoredMapList = ExSequencing(M,N,ants,acoIteration,finalLSIteration, 3, testCounter)
            subString = ProcessResult(highestScoredMapList, realPeptide)
            resultString = resultString + subString
            print("subString = " + subString)
            print(resultString)
            highestScoredMapList = ExSequencing(M,N,ants,acoIteration,finalLSIteration, 4, testCounter)
            subString = ProcessResult(highestScoredMapList, realPeptide)
            resultString = resultString + subString
            resultString = resultString + "\n"
            print("subString = " + subString)
            print(resultString)

        resultString = resultString + "\n"
    print("resultString = ")
    print(resultString)
#Batch test segment ends here

#Scoring segment starts here
def GetBestScoredPeptidesFromMap(comparisonMap,peptideScore,considerEqual):
    keys = list(comparisonMap.keys())
    keys.sort(reverse = True)
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


def SequenceScoreExt(peptide):
    score = 0
    sequenceBreakDownMap = {}
    subsequenceMap = {}
    SequenceSpectrum(peptide, sequenceBreakDownMap, subsequenceMap)
    tmpSpectrumMap = copy(SpectrumMap)
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




def SequenceSpectrumExt(Peptide,sequenceBreakDownMap,subsequenceMap):
    PrefixMass = []
    Peptides = []
    Peptides = Peptide.split("-")
    Peptides = [int(x) for x in Peptides]
    n = len(Peptides)
    Peptides.extend(Peptides)
    for i in range(n):
        tmp = Peptides[i]
        tmpList = []
        tmpList.append(Peptides[i])
        if i == 0:
            for j in range(i+1,n):
                tmp += Peptides[j]
                tmpList.append(Peptides[j])
                if tmp not in sequenceBreakDownMap:
                    sequenceBreakDownMap[tmp] = []
                sequenceBreakDownMap[tmp].extend(tmpList)
        else:
            for j in range(i+1,n+i-1):
                tmp += Peptides[j]
                tmpList.append(Peptides[j])
                if tmp not in sequenceBreakDownMap:
                    sequenceBreakDownMap[tmp] = []
                sequenceBreakDownMap[tmp].extend(tmpList)

    for i in range(n):
        tmp = Peptides[i]
        if tmp not in subsequenceMap:
            subsequenceMap[tmp] = []
        tmpList = []
        if i == 0:
            for j in range(i+1,n):
                tmp += Peptides[j]
                tmpList.append(tmp)
            subsequenceMap[Peptides[i]].extend(tmpList)
        else:
            for j in range(i+1,n+i-1):
                tmp += Peptides[j]
                tmpList.append(tmp)
            subsequenceMap[Peptides[i]].extend(tmpList)

    if 0 not in subsequenceMap:
        subsequenceMap[0] = []


def SequenceSpectrum(Peptide,sequenceBreakDownMap,subsequenceMap):
    PrefixMass = []
    Peptides = []
    Peptides = Peptide.split("-")
    Peptides = [int(x) for x in Peptides]
    n = len(Peptides)
    for i in range(n):
        tmp = Peptides[i]
        tmpList = []
        tmpList.append(Peptides[i])
        for j in range(i+1,n):
            tmp += Peptides[j]
            tmpList.append(Peptides[j])
            if tmp not in sequenceBreakDownMap:
                sequenceBreakDownMap[tmp] = []
            sequenceBreakDownMap[tmp].extend(tmpList)
    for i in range(n):
        tmp = Peptides[i]
        if tmp not in subsequenceMap:
            subsequenceMap[tmp] = []
        tmpList = []
        for j in range(i+1,n):
            tmp += Peptides[j]
            tmpList.append(tmp)
        subsequenceMap[Peptides[i]].extend(tmpList)
    if 0 not in subsequenceMap:
        subsequenceMap[0] = []



def SequenceScore(peptide):
    score = 0
    sequenceBreakDownMap = {}
    subsequenceMap = {}
    SequenceSpectrum(peptide, sequenceBreakDownMap, subsequenceMap)
    keys = list(subsequenceMap.keys())
    for i in range(len(keys)):
        subsequence = keys[i]
        subsequenceList = subsequenceMap[keys[i]]
        if subsequence in SpectrumMap:
            score += 1
            for j in range(len(subsequenceList)):
                if subsequenceList[j] in SpectrumMap:
                    score += 2

    return score


def LinearScore(peptide):
    score = 0
    testSpectrum = LinearSpectrum(peptide)
    prevSubpeptide = ""
    for subpeptide in testSpectrum:
        if subpeptide == prevSubpeptide:
            continue
        subpeptideCount = testSpectrum.count(subpeptide)
        try:
            tmp = SpectrumMap[subpeptide]
        except Exception:
            tmp = 0
        if  tmp <= subpeptideCount and tmp > 0:
            score += tmp
        elif tmp > subpeptideCount:
            score += subpeptideCount
        prevSubpeptide = subpeptide
    return score


def CycloScore(peptide):
    score = 0
    if peptide == "":
        return 0
    testSpectrum = CyclicSpectrum(peptide)
    prevSubpeptide = ""

    for subpeptide in testSpectrum:
        if subpeptide == prevSubpeptide:
            continue
        subpeptideCount = testSpectrum.count(subpeptide)
        try:
            tmp = SpectrumMap[subpeptide]
        except Exception:
            tmp = 0
        if  tmp <= subpeptideCount and tmp > 0:
            score += tmp
        elif tmp > subpeptideCount:
            score += subpeptideCount
        prevSubpeptide = subpeptide
    return score


def LinearSpectrum(Peptide):
    linearSpectrum = []
    PrefixMass = []
    Peptides = []
    if type(Peptide) is np.int32:
        Peptides.append(Peptide)
        n = 1
    else:
        Peptides = Peptide.split("-")
        n = len(Peptides)
    for i in range(n+1):
        PrefixMass.append(0)
    for i in range(n):
        mass = int(Peptides[i])
        PrefixMass[i+1] = PrefixMass[i] + mass
    for i in range(n):
        for j in range(i+1,n+1):
            tmp = PrefixMass[j] - PrefixMass[i]
            linearSpectrum.append(tmp)
    linearSpectrum.append(0)
    linearSpectrum.sort()
    return linearSpectrum


def CyclicSpectrum(Peptide):
    cyclicSpectrum = []
    PrefixMass = []
    Peptides = []
    if type(Peptide) is np.int32:
        Peptides.append(Peptide)
        n = 1
    else:
        Peptides = Peptide.split("-")
        n = len(Peptides)
    for i in range(n+1):
        PrefixMass.append(0)
    for i in range(n):
        mass = int(Peptides[i])
        PrefixMass[i+1] = PrefixMass[i] + mass
    peptideMass = PrefixMass[n]
    for i in range(n):
        for j in range(i+1,n+1):
            tmp = PrefixMass[j] - PrefixMass[i]
            cyclicSpectrum.append(tmp)
            if i > 0 and j < n :
                tmp2 = peptideMass - tmp
                cyclicSpectrum.append(tmp2)
    cyclicSpectrum.append(0)
    cyclicSpectrum.sort()
    return cyclicSpectrum

#Scoring segment ends here

#operator segment starts
def GetBestScoredPeptidesFromMap(comparisonMap,peptideScore,considerEqual):
    keys = list(comparisonMap.keys())
    keys.sort(reverse = True)
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


def CalculateMassOfPeptide(peptide):
    peptide = peptide.strip("-")
    peptideList = peptide.split('-')
    peptideList = [int(i) for i in peptideList]
    return sum(peptideList)


def Fusion(peptide,extensiveScoring,fitnessObj):
    global FussionSuccess
    splittedPeptide = copy(peptide.peptide.split("-"))
    splittedPeptide = [int(x) for x in splittedPeptide]
    concatKeys = list(MapToConcatenate.keys())
    isChanged = False
    resultlist = []
    comparisonMap = {}

    for firstFragmentIndex in range(len(splittedPeptide)):
        if splittedPeptide[firstFragmentIndex] in concatKeys:
            concatValues = MapToConcatenate[splittedPeptide[firstFragmentIndex]]
            for secondFragment in range(len(concatValues)):
                if concatValues[secondFragment] in splittedPeptide:
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
                        generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
                    else:
                        generatedPeptide.score = SequenceScore(generatedPeptide.peptide)
                    fitnessObj.increaseFE(generatedPeptide)

                    if generatedPeptide.score > peptide.score:
                        FussionSuccess += 1
                        if generatedPeptide.score not in comparisonMap:
                            comparisonMap[generatedPeptide.score] = []
                        comparisonMap[generatedPeptide.score].append(generatedPeptide)
    bestScoredPeptideList = []
    if len(comparisonMap) > 0:
        bestScoredPeptideList = GetBestScoredPeptidesFromMap(comparisonMap,peptide.score,False)
        for k in range(len(bestScoredPeptideList)):
            bestScoredPeptideList[k].mass = CalculateMassOfPeptide(bestScoredPeptideList[k].peptide)

    return bestScoredPeptideList


def Fission(peptide,extensiveScoring,fitnessObj):
    global FissonSuccess
    splittedPeptide = copy(peptide.peptide.split("-"))
    splittedPeptide = [int(x) for x in splittedPeptide]
    splitKeys = list(MapToSplit.keys())
    isChanged = False
    resultlist = []
    comparisonMap = {}

    for i in range(len(splittedPeptide)):
        if splittedPeptide[i] in splitKeys:
            splitValues = MapToSplit[splittedPeptide[i]]
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
                    generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
                else:
                    generatedPeptide.score = SequenceScore(generatedPeptide.peptide)
                fitnessObj.increaseFE(generatedPeptide)

                if generatedPeptide.score > peptide.score:
                    FissonSuccess += 1
                    if generatedPeptide.score not in comparisonMap:
                        comparisonMap[generatedPeptide.score] = []
                    comparisonMap[generatedPeptide.score].append(generatedPeptide)
    bestScoredPeptideList = []
    if len(comparisonMap) > 0:
        bestScoredPeptideList = GetBestScoredPeptidesFromMap(comparisonMap,peptide.score,False)
        for k in range(len(bestScoredPeptideList)):
            bestScoredPeptideList[k].mass = CalculateMassOfPeptide(bestScoredPeptideList[k].peptide)
            bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)

    return bestScoredPeptideList


def Swap(peptide, extensiveScoring, fitnessObj):
    global SwapSuccess
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

        generatedPeptide = PeptideItem()
        for k in range(len(splittedPeptide)):
            generatedPeptide.peptide += str(splittedPeptide[k]) + "-"
        generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
        if extensiveScoring:
            generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
        else:
            generatedPeptide.score = SequenceScore(generatedPeptide.peptide)

        fitnessObj.increaseFE(generatedPeptide)
        if generatedPeptide.score > peptide.score:
            SwapSuccess += 1
            if generatedPeptide.score not in comparisonMap:
                comparisonMap[generatedPeptide.score] = []
            comparisonMap[generatedPeptide.score].append(generatedPeptide)

    bestScoredPeptideList = []
    if len(comparisonMap) > 0:
        bestScoredPeptideList = GetBestScoredPeptidesFromMap(comparisonMap,peptide.score,True)
        for k in range(len(bestScoredPeptideList)):
            bestScoredPeptideList[k].mass = CalculateMassOfPeptide(bestScoredPeptideList[k].peptide)
            bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)

    return bestScoredPeptideList


def Transposition(peptide,extensiveScoring,fitnessObj):
    global TransportationSuccess
    comparisonMap = {}

    for i in range(5):
        splittedPeptide = copy(peptide.peptide.split("-"))
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
            generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
        else:
            generatedPeptide.score = SequenceScore(generatedPeptide.peptide)
        fitnessObj.increaseFE(generatedPeptide)
        if generatedPeptide.score > peptide.score:
            TransportationSuccess += 1
            if generatedPeptide.score not in comparisonMap:
                comparisonMap[generatedPeptide.score] = []
            comparisonMap[generatedPeptide.score].append(generatedPeptide)

    bestScoredPeptideList = []
    if len(comparisonMap) > 0:
        bestScoredPeptideList = GetBestScoredPeptidesFromMap(comparisonMap,peptide.score,True)
        for k in range(len(bestScoredPeptideList)):
            bestScoredPeptideList[k].mass = CalculateMassOfPeptide(bestScoredPeptideList[k].peptide)
            bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)

    return bestScoredPeptideList



def Mutate(peptide,extensiveScoring, fitnessObj):
    subsequenceMap = {}
    sequenceBreakDownMap = {}
    SequenceSpectrum(peptide.peptide, sequenceBreakDownMap,subsequenceMap)
    comparisonMap = {}
    comparisonMap[peptide.score] = []
    comparisonMap[peptide.score].append(peptide)
    aminoAcidToChange = -1
    keys = list(sequenceBreakDownMap.keys())
    keys.sort(reverse = True)
    shouldStop = False
    for i in range(len(keys)):
        isKeyPresent = False
        if keys[i] in SpectrumMap:
            isKeyPresent = True
        subsequenceList = sequenceBreakDownMap[keys[i]]
        for j in range(len(subsequenceList)):
            if subsequenceList[j] not in SpectrumMap:
                sequenceList = subsequenceMap[subsequenceList[j]]
                noOfSequencesPresent = 0
                for k in range(len(sequenceList)):
                    if sequenceList[k] in SpectrumMap:
                        noOfSequencesPresent += 1
                if noOfSequencesPresent < int(len(sequenceList)*2/3):
                    aminoAcidToChange = subsequenceList[j]
                    shouldStop = True
                    break
        if shouldStop:
            break
    splittedPeptide = copy(peptide.peptide.split("-"))
    splittedPeptide = [int(x) for x in splittedPeptide]
    if aminoAcidToChange != -1:
        positionTochange = splittedPeptide.index(aminoAcidToChange)
    else:
        positionTochange = np.random.randint(0,len(splittedPeptide))

    for j in range(3):
        rand1 = np.random.randint(0,len(AminoAcidMass))
        splittedPeptide[positionTochange] = AminoAcidMass[rand1]
        generatedPeptide = PeptideItem()
        for k in range(len(splittedPeptide)):
            generatedPeptide.peptide += str(splittedPeptide[k]) + "-"
        generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
        if extensiveScoring:
            generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
        else:
            generatedPeptide.score = SequenceScore(generatedPeptide.peptide)

        fitnessObj.increaseFE(generatedPeptide)
        if generatedPeptide.score not in comparisonMap:
            comparisonMap[generatedPeptide.score] = []
        comparisonMap[generatedPeptide.score].append(generatedPeptide)
    bestScoredPeptideList = GetBestScoredPeptidesFromMap(comparisonMap,peptide.score,True)
    for k in range(len(bestScoredPeptideList)):
        bestScoredPeptideList[k].mass = CalculateMassOfPeptide(bestScoredPeptideList[k].peptide)
        bestScoredPeptideList[k].ringLength = len(bestScoredPeptideList[k].peptide)

    return bestScoredPeptideList



def Reverse(peptideItemList,extensiveScoring):
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
            generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
        else:
            generatedPeptide.score = SequenceScore(generatedPeptide.peptide)
        generatedPeptide.mass = CalculateMassOfPeptide(generatedPeptide.peptide)
        reversePeptideList.append(generatedPeptide)

    return reversePeptideList



def MergeAndBreak(peptide,extensiveScoring,fitnessObj):
    global MergeAndBreakSuccess
    comparisonMap = {}
    splittedPeptide = copy(peptide.peptide.split("-"))
    splittedPeptide = [int(x) for x in splittedPeptide]
    mapKeys = list(MergeAndBreakMap.keys())
    shouldBreak = False
    for i in range(len(splittedPeptide)-1):
        sumOfMass = splittedPeptide[i] + splittedPeptide[i+1]
        if sumOfMass in mapKeys:
            mapValues = MergeAndBreakMap[sumOfMass]
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
                        generatedPeptide.score = SequenceScoreExt(generatedPeptide.peptide)
                    else:
                        generatedPeptide.score = SequenceScore(generatedPeptide.peptide)
                    fitnessObj.increaseFE(generatedPeptide)
                    if generatedPeptide.score > peptide.score:
                        generatedPeptide.mass = CalculateMassOfPeptide(generatedPeptide.peptide)
                        shouldBreak = True
                        MergeAndBreakSuccess += 1
                        break
        if shouldBreak:
            break
    bestScoredPeptideList = []
    if shouldBreak:
        bestScoredPeptideList.append(generatedPeptide)
    return bestScoredPeptideList


#operator segments ends

#main segment starts
def CreateAminoAcidMergeAndBreakMap():
    global MergeAndBreakMap
    MergeAndBreakMap.clear()
    for i in range(len(AminoAcidMass)-1):
        firstAminoAcid = AminoAcidMass[i]
        for j in range(i+1,len(AminoAcidMass)):
            secondAminoAcid = AminoAcidMass[j]
            sumOfMass = firstAminoAcid + secondAminoAcid
            if sumOfMass not in MergeAndBreakMap:
                MergeAndBreakMap[sumOfMass] = []
            MergeAndBreakMap[sumOfMass].append([firstAminoAcid,secondAminoAcid])


def CreateAminoAcidRearrangementMap():
    global MapToConcatenate
    global MapToSplit
    MapToConcatenate.clear()
    MapToSplit.clear()

    for i in range(len(AminoAcidMass)):
        k = AminoAcidMass[i]
        for j in range(i+1,len(AminoAcidMass)):
            v = AminoAcidMass[j]
            sumOfMass = k + v
            if sumOfMass > 200:
                break
            if sumOfMass in AminoAcidMass:
                if k not in MapToConcatenate:
                    MapToConcatenate[k] = []
                MapToConcatenate[k].append(v)
                if v not in MapToConcatenate:
                    MapToConcatenate[v] = []
                MapToConcatenate[v].append(k)
                if sumOfMass not in MapToSplit:
                    MapToSplit[sumOfMass] = []
                MapToSplit[sumOfMass].append([k,v])


def GetFrequentAminoAcids(M):
    global AminoAcidMass
    AminoAcidMass.clear()
    output = []
    countMap ={}
    for i in range(len(Spectrum)):
        for j in range(i):
            diff = abs(Spectrum[j] - Spectrum[i])
            if diff >= 57 and diff <=200:
                output.append(diff)
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

    keys = list(countMap.keys())
    keys.sort(reverse = True)
    for index in range(len(keys)):
        tmpArr = countMap[keys[index]]
        AminoAcidMass.extend(tmpArr)
        if len(AminoAcidMass) >= M:
            break

    AminoAcidMass.sort()
    return AminoAcidMass


def PopulateSpectrumMap():
    global SpectrumMap
    global Spectrum
    SpectrumMap.clear()
    for i in range(len(Spectrum)):
        item = Spectrum[i]
        if item not in SpectrumMap:
            SpectrumMap[item] = 0
        SpectrumMap[item] += 1


def AcoLocal(ants,N, acoIteration,finalLSIteration):
    global ParentMass
    global Pheromone
    global AcoLocalFitnessEvaluation
    AcoLocalFitnessEvaluation = 0
    fitnessObj = FitnessEvaluation(1)

    print("in AcoLocal")
    print("acoIteration = " + str(acoIteration))
    print("finalLSIteration = " + str(finalLSIteration))
    HighestScoredPeptideMap.clear()
    noOfColumns = len(AminoAcidMass)
    noOfRows = np.ceil(ParentMass / AminoAcidMass[0])
    Pheromone = np.ones((noOfRows,noOfColumns),dtype = float)
    pheromoneDimension = Pheromone.shape
    pheromoneRows = pheromoneDimension[0]
    pheromoneColumns = pheromoneDimension[1]
    #aco
    for epoch in range(2):
        Pheromone.fill(1.)
        probability = Pheromone.copy()
        probability = probability + 0.1
        probability = probability/probability.sum(axis=1)[:,None]
        lenAminoAcidMass = len(AminoAcidMass)
        resultMap = {}
        for k in range(acoIteration):
            peptideItemList = []
            scoreList = []
            #print("no of iteration = " + str(k))
            for i in range(ants):
                newPeptide = CreateAntPath(probability,False)
                peptideItemList.append(newPeptide)
                fitnessObj.increaseFE(newPeptide)

            modifiedPeptideItemList = PerformLocalSearch(peptideItemList,fitnessObj)
            CalculatePheromone(modifiedPeptideItemList)
            SaveHighestScoredPeptides(modifiedPeptideItemList, N, False)
            probability = Pheromone.copy()
            probability = probability + 0.1
            probability = probability/probability.sum(axis=1)[:,None]

    print("HighestScoredPeptideMap after first 4000 iterations")
    for k,v in HighestScoredPeptideMap.items():
        print(str(k) + "-> [",end="")
        for j in v:
            print(j.peptide + "(" + str(j.mass) + "->" + str(j.score) + ")",end="")
            print(",")
        print("]")

    peptideItemList = GetHighestScoredPeptides(True)
    print("highestScored peptides after conversion to SequenceScoreExt")
    for i in range(len(peptideItemList)):
        print(peptideItemList[i].peptide + "-> (" + str(peptideItemList[i].score) + ")" )
    HighestScoredPeptideMap.clear()
    peptideItemList = ProcessPathsFurther(peptideItemList,finalLSIteration,N,False,1, fitnessObj)
    AcoLocalFitnessEvaluation = fitnessObj.getFitnessEvaluationCount()
    fitnessObj.printFE()
    print("AcoLocalFitnessEvaluation = " + str(AcoLocalFitnessEvaluation))
    print("FitnessEvaluation = " + str(fitnessObj.getFitnessEvaluationCount()))
    print("HighestScoredPeptideMap after last 200 iterations")
    for k,v in HighestScoredPeptideMap.items():
        print(str(k) + "-> [",end="")
        for j in v:
            print(j.peptide + "(" + str(j.mass) + "->" + str(j.score) + ")",end="")
            print(",")
        print("]")


def ProcessPathsFurther(peptideItemList, counter,N, isLocalSearchOnly, method, fitnessObj):
    global AcoLocalFitnessEvaluation
    modifiedPeptideItemList = copy(peptideItemList)
    j = 0
    if method > 1:
        j = fitnessObj.getFitnessEvaluationCount()
        counter = AcoLocalFitnessEvaluation
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
                resultlist = Mutate(tmpPeptide,True,fitnessObj)
            elif randomOperator == 2:
                resultlist = Transposition(tmpPeptide,True,fitnessObj)
            elif randomOperator == 3:
                resultlist = Swap(tmpPeptide,True,fitnessObj)
            elif randomOperator == 4:
                resultlist = Fusion(tmpPeptide,True,fitnessObj)
            elif randomOperator == 5:
                resultlist = Fission(tmpPeptide,True,fitnessObj)
            else:
                resultlist = MergeAndBreak(tmpPeptide,True,fitnessObj)
            if len(resultlist) > 0:
                modifiedPeptideItemList.extend(resultlist)
            i += 1

        SaveHighestScoredPeptides(modifiedPeptideItemList, N, isLocalSearchOnly)
        modifiedPeptideItemList = GetHighestScoredPeptides(False)
        '''print("modifiedPeptideItemList after path process j = " + str(j))
        for k in modifiedPeptideItemList:
            print(k.peptide + "(" + str(k.score) + " -> "+str(k.mass)+")")'''
        if method > 1:
            j = fitnessObj.getFitnessEvaluationCount()
        else:
            j += 1
    reversePeptideList = Reverse(modifiedPeptideItemList,True)
    modifiedPeptideItemList.extend(reversePeptideList)
    SaveHighestScoredPeptides(modifiedPeptideItemList, N, isLocalSearchOnly)

    return modifiedPeptideItemList



def GetHighestScoredPeptides(withNewScore):
    highestScoredPeptideList = []
    keys = HighestScoredPeptideMap.keys()
    for i in keys:
        values = HighestScoredPeptideMap[i]
        for peptide in values:
            if withNewScore:
                peptide.score = SequenceScoreExt(peptide.peptide)
            highestScoredPeptideList.append(peptide)

    return highestScoredPeptideList



def SaveHighestScoredPeptides(peptideItemList,N, isLocalSearchOnly):
    global ParentMass
    peptideScoreMap = {}
    for i in range(len(peptideItemList)):
        tmpPeptide = peptideItemList[i]
        if isLocalSearchOnly:
            if tmpPeptide.mass in range(ParentMass-30, ParentMass+30):
                if tmpPeptide.score not in peptideScoreMap:
                    peptideScoreMap[tmpPeptide.score] = []
                peptideScoreMap[tmpPeptide.score].append(tmpPeptide)
        else:
            if tmpPeptide.mass == ParentMass:
                if tmpPeptide.score not in peptideScoreMap:
                    peptideScoreMap[tmpPeptide.score] = []
                peptideScoreMap[tmpPeptide.score].append(tmpPeptide)

    mergedKeys = list(HighestScoredPeptideMap.keys())
    peptideScoreMapKeys = list(peptideScoreMap.keys())
    mergedKeys.extend(peptideScoreMapKeys)
    mergedKeys = np.unique(mergedKeys)
    mergedKeys[::-1].sort()
    count = 0
    tmpMap = copy(HighestScoredPeptideMap)
    HighestScoredPeptideMap.clear()
    for i in range(len(mergedKeys)):
        k = mergedKeys[i]
        v = tmpMap.get(k,[])
        v.extend(peptideScoreMap.get(k,[]))
        HighestScoredPeptideMap[k] = []
        values = []
        for j in range(len(v)):
            if v[j].peptide not in values:
                values.append(v[j].peptide)
                HighestScoredPeptideMap[k].append(v[j])
        count += len(values)
        if count >= N:
            break


def PerformLocalSearch(peptideItemList,fitnessObj):
    global LocalSearchIteration
    '''print("incoming peptideItemList")
    for i in peptideItemList:
        print(i.peptide + "(" + str(i.score) + " -> "+str(i.mass)+")")'''
    modifiedPeptideItemList = copy(peptideItemList)
    for m in range(LocalSearchIteration):
        i = 0
        while i < len(modifiedPeptideItemList):
            tmpPeptide = modifiedPeptideItemList[i]
            resultlist = Mutate(tmpPeptide, False,fitnessObj)
            if len(resultlist) > 0:
                del modifiedPeptideItemList[i]
                for k in range(len(resultlist)):
                    modifiedPeptideItemList.insert(i+k,resultlist[k])
                i += len(resultlist)
            else:
                i += 1

    '''print("outgoing modifiedPeptideItemList")
    for k in modifiedPeptideItemList:
        print(k.peptide + "(" + str(k.score) + " -> "+str(k.mass)+")")'''
    return modifiedPeptideItemList


def CreateAntPath(probability,extensiveScoring):
    global ParentMass
    newPeptide = PeptideItem()

    while newPeptide.mass <= ParentMass:
        if (ParentMass-newPeptide.mass) < AminoAcidMass[0]:
            break
        splittedNewPeptide = newPeptide.peptide.split("-")
        peptideLen = len(splittedNewPeptide)
        probabilityRow = probability[peptideLen-1].copy()
        dimensionProbability = probability.shape
        massToBeAdded = np.random.choice(AminoAcidMass,p=probabilityRow)
        newPeptide.peptide += str(massToBeAdded) + "-"
        newPeptide.mass = CalculateMassOfPeptide(newPeptide.peptide)
    #remove last -
    newPeptide.peptide = newPeptide.peptide.strip("-")
    if extensiveScoring:
        newPeptide.score = SequenceScoreExt(newPeptide.peptide)
    else:
        newPeptide.score = SequenceScore(newPeptide.peptide)
    return newPeptide




def CalculatePheromone(peptideItemList):

    global Pheromone
    evaporateQuantity = 1 - PheromoneConstant
    Pheromone = Pheromone * evaporateQuantity
    # add pheromone for all ants
    for i in range(len(peptideItemList)):
        if peptideItemList[i].score > 0:
            tmpPeptide = peptideItemList[i].peptide
            aminoAcidsInTmpPeptide = tmpPeptide.split('-')
            aminoAcidsInTmpPeptide = [int(j) for j in aminoAcidsInTmpPeptide]
            pheromoneToBeAdded = peptideItemList[i].score * Q
            for j in range(len(aminoAcidsInTmpPeptide)):
                columnIndexPheromonMatrix = AminoAcidMass.index(aminoAcidsInTmpPeptide[j])
                rowIndexInPheromoneMatrix = j
                Pheromone[rowIndexInPheromoneMatrix,columnIndexPheromonMatrix] += pheromoneToBeAdded


def AcoOnly(ants,N, acoIteration,finalLSIteration):
    global ParentMass
    global Pheromone
    global AcoLocalFitnessEvaluation
    fitnessObj = FitnessEvaluation(2)
    print("acoIteration = " + str(acoIteration))
    print("finalLSIteration = " + str(finalLSIteration))
    print("AcoOnly")
    HighestScoredPeptideMap.clear()
    noOfColumns = len(AminoAcidMass)
    noOfRows = np.ceil(ParentMass / AminoAcidMass[0])

    Pheromone = np.ones((noOfRows,noOfColumns),dtype = float)
    pheromoneDimension = Pheromone.shape
    pheromoneRows = pheromoneDimension[0]
    pheromoneColumns = pheromoneDimension[1]
    #aco
    half = int(AcoLocalFitnessEvaluation / 2)
    print("half in AcoOnly = " + str(half))
    for epoch in range(2):
        Pheromone.fill(1.)
        probability = Pheromone.copy()
        probability = probability + 0.1
        probability = probability/probability.sum(axis=1)[:,None]
        lenAminoAcidMass = len(AminoAcidMass)
        extensiveScoring = False
        #2050
        while fitnessObj.getFitnessEvaluationCount() < AcoLocalFitnessEvaluation:
            peptideItemList = []
            #print("no of iteration = " + str(k))
            for i in range(ants):
                newPeptide = CreateAntPath(probability,False)
                peptideItemList.append(newPeptide)
                fitnessObj.increaseFE(newPeptide)

            CalculatePheromone(peptideItemList)
            SaveHighestScoredPeptides(peptideItemList,N, extensiveScoring)
            probability = Pheromone.copy()
            probability = probability + 0.1
            probability = probability/probability.sum(axis=1)[:,None]
            if fitnessObj.getFitnessEvaluationCount() == half:
                extensiveScoring = True
                break

    fitnessObj.printFE()
    print("AcoLocalFitnessEvaluation = " + str(AcoLocalFitnessEvaluation))
    print("FitnessEvaluation = " + str(fitnessObj.getFitnessEvaluationCount()))

    print("HighestScoredPeptideMap after first 2000 iterations")
    for k,v in HighestScoredPeptideMap.items():
        print(str(k) + "-> [",end="")
        for j in v:
            print(j.peptide + "(" + str(j.mass) + "->" + str(j.score) + ")",end="")
            print(",")
        print("]")


def LocalSearchOnly(N,ants,acoIteration,finalLSIteration):
    global ParentMass
    global AcoLocalFitnessEvaluation
    fitnessObj = FitnessEvaluation(3)
    HighestScoredPeptideMap.clear()
    print("acoIteration = " + str(acoIteration))
    print("finalLSIteration = " + str(finalLSIteration))
    print("LocalSearchOnly")

    randomlyGeneratedPeptide = []
    for i in range(ants):
        generatedPeptide = PeptideItem()
        while generatedPeptide.mass <= ParentMass:
            rand = np.random.randint(0,len(AminoAcidMass))
            generatedPeptide.peptide += str(AminoAcidMass[rand]) + "-"
            generatedPeptide.mass += AminoAcidMass[rand]
        generatedPeptide.peptide = generatedPeptide.peptide.strip("-")
        generatedPeptide.score = SequenceScore(generatedPeptide.peptide)
        randomlyGeneratedPeptide.append(generatedPeptide)
        fitnessObj.increaseFE(generatedPeptide)

    peptideItemList = copy(randomlyGeneratedPeptide)
    for counter in range(2):
        for k in range(acoIteration):
            if len(peptideItemList) > 0:
                modifiedPeptideItemList = PerformLocalSearch(peptideItemList,fitnessObj)
                SaveHighestScoredPeptides(modifiedPeptideItemList,N, True)
                peptideItemList = GetHighestScoredPeptides(False)

            else:
                peptideItemList = randomlyGeneratedPeptide
                break
    peptideItemList = GetHighestScoredPeptides(True)
    HighestScoredPeptideMap.clear()
    peptideItemList = ProcessPathsFurther(peptideItemList,finalLSIteration,N, True, 3, fitnessObj)

    fitnessObj.printFE()
    print("AcoLocalFitnessEvaluation = " + str(AcoLocalFitnessEvaluation))
    print("FitnessEvaluation = " + str(fitnessObj.getFitnessEvaluationCount()))

    print("HighestScoredPeptideMap after last 50 iterations")
    for k,v in HighestScoredPeptideMap.items():
        print(str(k) + "-> [",end="")
        for j in v:
            print(j.peptide + "(" + str(j.mass) + "->" + str(j.score) + ")",end="")
            print(",")
        print("]")



def BnB(N, acoIteration):
    global ParentMass
    HighestScoredPeptideMap.clear()
    global AcoLocalFitnessEvaluation
    fitnessObj = FitnessEvaluation(4)
    LeaderBoard = []
    LeaderBoard.append("")

    LeaderPeptide = ""
    MaxScore = 0
    counter = 0
    maxScoredPeptide = []
    generatedPeptide = PeptideItem()

    while fitnessObj.getFitnessEvaluationCount() < AcoLocalFitnessEvaluation and len(LeaderBoard) > 0 :
        LeaderBoard = Expand(LeaderBoard)
        tmpLeaderBoard =  deepcopy(LeaderBoard)
        for i in range(len(LeaderBoard)):
            peptide = LeaderBoard[i]
            mass = Mass(peptide)
            if mass == ParentMass:
                cyscore = CycloScore(peptide)
                LeaderPeptideCyScore = CycloScore(LeaderPeptide)
                generatedPeptide.peptide = peptide
                generatedPeptide.score = cyscore
                fitnessObj.increaseFE(generatedPeptide)
                generatedPeptide.peptide = LeaderPeptide
                generatedPeptide.score = LeaderPeptideCyScore
                fitnessObj.increaseFE(generatedPeptide)

                if cyscore >= LeaderPeptideCyScore:
                    LeaderPeptide = peptide
                    if cyscore > MaxScore:
                        maxScoredPeptide = []
                        maxScoredPeptide.append(peptide)
                    elif cyscore == MaxScore:
                        maxScoredPeptide.append(peptide)
                    MaxScore = cyscore
            elif mass > ParentMass:
                tmpLeaderBoard.remove(peptide)

        if len(tmpLeaderBoard) > 0:
            LeaderBoard = Trim(tmpLeaderBoard,N, fitnessObj)

        else:
            LeaderBoard = tmpLeaderBoard
        counter += 1
    print("counter = " + str(counter))

    modifiedPeptideItemList = []
    for i in range(len(maxScoredPeptide)):
        score = SequenceScoreExt(maxScoredPeptide[i])
        generatedPeptide = PeptideItem()
        generatedPeptide.peptide = maxScoredPeptide[i]
        generatedPeptide.score = score
        generatedPeptide.mass = CalculateMassOfPeptide(maxScoredPeptide[i])
        modifiedPeptideItemList.append(generatedPeptide)
    SaveHighestScoredPeptides(modifiedPeptideItemList, N, False)

    fitnessObj.printFE()
    print("AcoLocalFitnessEvaluation = " + str(AcoLocalFitnessEvaluation))
    print("FitnessEvaluation = " + str(fitnessObj.getFitnessEvaluationCount()))

    print("HighestScoredPeptideMap in BnB")
    for k,v in HighestScoredPeptideMap.items():
        print(str(k) + "-> [",end="")
        for j in v:
            print(j.peptide + "(" + str(j.mass) + "->" + str(j.score) + ")",end="")
            print(",")
        print("]")




def Mass(peptide):
    mass = 0
    try:
        peptides = peptide.split('-')
        for i in range(len(peptides)):
            mass += int(peptides[i])
    except Exception:
        mass = peptide
    return mass


def Expand(Peptides):
    extendedPeptide = []
    if Peptides[0] == "":
        extendedPeptide = AminoAcidMass
    else:
        for i in range(len(Peptides)):
            for j in range(len(AminoAcidMass)):
                tmpPeptide = str(Peptides[i]) + '-' +str(AminoAcidMass[j])
                extendedPeptide.append(tmpPeptide)
    return extendedPeptide


def Trim(Leaderboard, N, fitnessObj):
    LeaderboardLen = len(Leaderboard)
    LinearScoreMap = {}
    LinearScoreList = []
    generatedPeptide = PeptideItem()
    for i in range(LeaderboardLen):
        peptide = Leaderboard[i]
        score = LinearScore(peptide)
        generatedPeptide.peptide = peptide
        generatedPeptide.score = score
        fitnessObj.increaseFE(generatedPeptide)
        LinearScoreList.append(score)
        if score not in LinearScoreMap:
            LinearScoreMap[score] = []
        LinearScoreMap[score].append(peptide)
    LinearScoreList.sort(reverse = True)
    for j in range(N,len(LinearScoreList)):
        if LinearScoreList[j] < LinearScoreList[N-1]:
            del LinearScoreList[(j):]
            break

    Leaderboard = []
    peptides = LinearScoreMap[LinearScoreList[0]]
    for j in range(len(peptides)):
        Leaderboard.append(peptides[j])
    for i in range(1,len(LinearScoreList)):
        if LinearScoreList[i] == LinearScoreList[i-1]:
            continue
        peptides = LinearScoreMap[LinearScoreList[i]]
        for j in range(len(peptides)):
            Leaderboard.append(peptides[j])
    return Leaderboard



def Sequencing(N,ants,acoIteration,finalLSIteration,method):
    global ParentMass
    global Spectrum
    ParentMass = max(Spectrum)
    #1 for aco local
    #2 for only aco
    #3 for only local search
    #4 for BnB
    if method == 1:
        AcoLocal(ants,N,acoIteration,finalLSIteration)
    elif method == 2:
        AcoOnly(ants,N,acoIteration,finalLSIteration)
    elif method == 3:
        LocalSearchOnly(ants,N,acoIteration,finalLSIteration)
    else:
        BnB(N,acoIteration)


def ExSequencing(M,N,ants,acoIteration,finalLSIteration, method, testCounter):
    resultList = []
    for i in range(testCounter):
        Sequencing(N,ants,acoIteration,finalLSIteration,method)
        resultMap = copy(HighestScoredPeptideMap)
        resultList.append(resultMap)
        fsock.flush()
    return resultList

def CalculateAvgFEScore(score):
    scoreAvgList = []
    for i in range(len(score[0])):
        tmpSum = 0
        for j in range(len(score)):
            tmpSum += score[j][i]
        tmpAvg = tmpSum / len(score)
        scoreAvgList.append(tmpAvg)
    return scoreAvgList


def main():
    """Main entry point for the script."""
    global Spectrum
    M,N,ants,Reads = sys.stdin.read().splitlines()
    M = int(M)
    N = int(N)
    ants = int(ants)
    Reads = Reads.split(",")
    Spectrum = [int(i) for i in Reads]

    print("Ant colony optimization output : ")
    print("Spectrum")
    print(Spectrum)
    print("M = " + str(M))
    print("N = " + str(N))
    print("ants = " + str(ants))
    PopulateSpectrumMap()
    GetFrequentAminoAcids(M)
    CreateAminoAcidRearrangementMap()
    CreateAminoAcidMergeAndBreakMap()
    Sequencing(N,ants,1600,200,1)
    fsock.flush()
    Sequencing(N,ants,1600,200,2)
    fsock.flush()
    Sequencing(N,ants,1600,200,3)
    fsock.flush()
    Sequencing(N,ants,1600,200,4)
    fsock.flush()


if __name__ == '__main__':
    fsock = open("E:\\bioInformatics\\aco_peptide\\output.txt", 'w')
    saveout = sys.stdout
    sys.stdout = fsock
    start = datetime.datetime.now()
    main()
    #Batchtest()
    end = datetime.datetime.now()
    print("total time required = " + str(end-start))
    sys.stdout = saveout
    fsock.close()
    sys.exit()
#main segment ends
