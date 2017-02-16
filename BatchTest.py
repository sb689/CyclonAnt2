from copy import copy,deepcopy
from Score import Score
from Init import Init
from PeptideSequence import PeptideSequence

class BatchTest:
    def __init__(self):
        self.scoreObj = None

    def setScoreObject(self,initObj):
        self.scoreObj = Score(initObj.spectrumMap)

    def selectParameters(self,spectrumLen):
        acoIteration = 1000
        finalLSIteration = 100
        M = 20
        N = 10
        ants = 10
        testCounter = 5
        if spectrumLen < 60:
            acoIteration = 1000
            finalLSIteration = 100
        elif spectrumLen in range(61,100):
            acoIteration = 1400
            finalLSIteration = 100
        elif spectrumLen in range(101,200):
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
        result.append(ants)
        result.append(testCounter)
        return result


    def getOriginalPeptides(self):
        originalPeptides = [0,0,0,0,0,0]
        originalPeptides[0] = "97-163-99-97-113-113-113-97"
        originalPeptides[1] = "99-114-113-147-97-99-114-113-147-97"
        originalPeptides[2] = "99-114-113-147-97-147-147-114-128-163"
        originalPeptides[3] = "99-114-113-147-97-163-99-114-113-147-97-163"
        originalPeptides[4] = "97-115-129-163-115-163-87-115-113-129-115-71-115"
        originalPeptides[5] = "99-114-113-147-97-163-99-114-113-147-97-163-88-96-174"
        return originalPeptides


    def getTestSpectrums(self,tag):
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



    def processResult(self,highestScoredMapList, realPeptide):
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
                    isMatch = matchPeptide(realPeptide, generatedPeptide.peptide)
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
            linearScore = self.scoreObj.linearScore(peptideToBeScored)
            cycloScore = self.scoreObj.cycloScore(peptideToBeScored)
            sequenceScore = self.scoreObj.sequenceScore(peptideToBeScored)
            sequenceScoreExt = self.scoreObj.sequenceScoreExt(peptideToBeScored)
            matched = 0
            if isMatch:
                matched = 1
            scoreList.append(cycloScore)
            returnString = returnString + str(linearScore) + ";" + str(cycloScore) + ";" + str(sequenceScore) + ";" + str(sequenceScoreExt) + ";" + str(matched) + ";"
        avgScore = sum(scoreList)/len(scoreList)
        maxScore = max(scoreList)
        returnString = str(maxScore) + ";" + str(avgScore) + ";" + returnString
        return returnString



def matchPeptide(originalPeptide, generatedPeptide):
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