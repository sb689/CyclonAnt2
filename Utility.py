class Utility:
    def __init__(self):
        self.method = 1

    def calculateMassOfPeptide(self,peptide):
        mass = 0
        peptide = peptide.strip("-")
        try:
            peptideList = peptide.split('-')
            peptideList = [int(i) for i in peptideList]
            mass = sum(peptideList)
        except Exception:
            mass = int(peptide)
        return mass
