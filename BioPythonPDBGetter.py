from Bio import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

class ResiGraph:
    Name = ""
    AdjList = {}
    AtomLocationsInAdjMatrix = {}
    def __init__(self, n, adjl):
        self.Name = n
        self.AtomLocationsInAdjMatrix = adjl

    def AddToAtomLocations(self, AtomName, Location):
        self.AddToAtomLocations[AtomName] = Location

    def AddToAdjList(self, Root, Neighbors):
        Neighbors = Neighbors
        self.AtomLocationsInAdjMatrix[Root, Neighbors]

ResiNum = 1
#Include allllll possible bonds just for graphing purposes

def ConstructGuanineConsecutive(GuanineResidue):
    ToReturn = {}
    GuanineConsecutiveAtoms = {"N9": ["C8", "C4", "C1'"],
                               "C8": ["H8", "N7", "N9"],
                               "H8": ["C8"],
                               "N7": ["C5", "C8"],
                               "C4": ["N3", "C5", "N9"],
                               "C6": ["O6", "N1", "C5"],
                               "O6": ["C6"],
                               "N1": ["H1", "C2", "C6"],
                               "H1": ["N1"],
                               "C2": ["N2", "N3", "N1"],
                               "N2": ["H21", "H22", "C2"],
                               "H21": ["N2"],
                               "H22": ["N2"],
                               "N3": ["C2", "C4"],
                               "C4": ["C5", "N3", "N9"],
                               "P": ["OP1", "O5'", "OP2"],
                               "OP1": ["P"],
                               "OP2": ["P"],
                               "O5'": ["P", "C5'"],
                               "C5'": ["H5'", "H5''", "C4'", "O5'"],
                               "C4'": ["H4'", "C3'", "O4'", "C5'"],
                               "C3'": ["O3'", "H3'", "C2'", "C4'"],
                               "O3'": ["C3'"],
                               "O4'": ["C4'", "C1'"],
                               "C2'": ["H2'", "C1'", "O2'", "C3'"],
                               "O2'": ["HO2'", "C2'"],
                               "H1'": ["C1'"],
                               "H2'": ["C2'"],
                               "H3'": ["C3'"],
                               "H4'": ["C4'"],
                               "H5'": ["C5'"],
                               "H5''": ["C5'"],
                               "HO2'": ["O2'"],
                               "C1'": ["O4'", "H1'", "C2'", "N9"]}

    for atom in GuanineResidue:
        AtomNameClean = atom.get_fullname().replace(" ", "").replace("''", "'")
        if AtomNameClean in GuanineConsecutiveAtoms:
            AtomsBondedToCurrentAtom = GuanineConsecutiveAtoms[AtomNameClean]
            VerticesToAdd = []
            for i in range(len(AtomsBondedToCurrentAtom)):
                try:
                    VerticesToAdd.append(GuanineResidue[AtomsBondedToCurrentAtom[i]])
                except:
                    pass
            ToReturn[atom] = VerticesToAdd
    return ToReturn

def ConstructAdenineConsecutive(AdenineResidue):
    ToReturn = {}
    AdenineConsecutiveAtoms = {"N9": ["C4", "C8", "C1'"],
                               "C8": ["N7", "H8", "N9"],
                               "C5": ["C6", "C4", "N7"],
                               "C6": ["N1", "N6", "C5"],
                               "N6": ["H61", "H62", "C6"],
                               "C2": ["H2", "N3", "N1"],
                               "N3": ["C4", "C2"],
                               "N1": ["C6", "C2"],
                               "C4": ["C5", "N3", "N9"],
                               "N7": ["C8", "C5"],
                               "H61": ["N6"],
                               "H62": ["N6"],
                               "H2": ["C2"],
                               "H8": ["C8"],
                               "P": ["OP1", "O5'", "OP2"],
                               "OP1": ["P"],
                               "OP2": ["P"],
                               "O5'": ["P", "C5'"],
                               "C5'": ["H5'", "H5''", "C4'", "O5'"],
                               "C4'": ["H4'", "C3'", "O4'", "C5'"],
                               "C3'": ["O3'", "H3'", "C2'", "C4'"],
                               "O3'": ["C3'"],
                               "O4'": ["C4'", "C1'"],
                               "C2'": ["H2'", "C1'", "O2'", "C3'"],
                               "O2'": ["HO2'", "C2'"],
                               "H1'": ["C1'"],
                               "H2'": ["C2'"],
                               "H3'": ["C3'"],
                               "H4'": ["C4'"],
                               "H5'": ["C5'"],
                               "H5''": ["C5'"],
                               "HO2'": ["O2'"],
                               "C1'": ["O4'", "H1'", "C2'", "N9"]}

    # Find out out HO2' is
    # Add hydrogens degree 1
    for atom in AdenineResidue:
        AtomNameClean = atom.get_fullname().replace(" ", "").replace("''", "'")
        if AtomNameClean in AdenineConsecutiveAtoms:
            AtomsBondedToCurrentAtom = AdenineConsecutiveAtoms[AtomNameClean]
            VerticesToAdd = []
            for i in range(len(AtomsBondedToCurrentAtom)):
                try:
                    VerticesToAdd.append(AdenineResidue[AtomsBondedToCurrentAtom[i]])
                except:
                    pass
            ToReturn[atom] = VerticesToAdd
    return ToReturn

def ConstructCytosineConsecutive(CytosineResidue):
    ToReturn = {}
    CytosineConsecutiveAtoms = {"N1": ["C2", "C6", "C1'"],
                               "C2": ["N3", "O2", "N1"],
                               "C4": ["N4", "C5", "N3"],
                               "N4": ["H41", "H42", "C4"],
                               "C5": ["H5", "C6", "C4"],
                                "N3": ["C2", "C4"],
                                "O2": ["C2"],
                                "C6": ["H6", "C5", "N1"],
                                "H41": ["N4"],
                                "H42": ["N4"],
                                "H5": ["C5"],
                                "H6": ["C6"],
                               "P": ["OP1", "O5'", "OP2"],
                               "OP1": ["P"],
                               "OP2": ["P"],
                               "O5'": ["P", "C5'"],
                               "C5'": ["H5'", "H5''", "C4'", "O5'"],
                               "C4'": ["H4'", "C3'", "O4'", "C5'"],
                               "C3'": ["O3'", "H3'", "C2'", "C4'"],
                               "O3'": ["C3'"],
                               "O4'": ["C4'", "C1'"],
                               "C2'": ["H2'", "C1'", "O2'", "C3'"],
                               "O2'": ["HO2'", "C2'"],
                               "H1'": ["C1'"],
                               "H2'": ["C2'"],
                               "H3'": ["C3'"],
                               "H4'": ["C4'"],
                               "H5'": ["C5'"],
                               "H5''": ["C5'"],
                               "HO2'": ["O2'"],
                               "C1'": ["O4'", "H1'", "C2'", "N1"]}


    # Find out out HO2' is
    # Add hydrogens degree 1
    for atom in CytosineResidue:
        AtomNameClean = atom.get_fullname().replace(" ", "").replace("''", "'")
        if AtomNameClean in CytosineConsecutiveAtoms:
            AtomsBondedToCurrentAtom = CytosineConsecutiveAtoms[AtomNameClean]
            VerticesToAdd = []
            for i in range(len(AtomsBondedToCurrentAtom)):
                try:
                    VerticesToAdd.append(CytosineResidue[AtomsBondedToCurrentAtom[i]])
                except:
                    pass
            ToReturn[atom] = VerticesToAdd
    return ToReturn

def ConstructUracilConsecutive(UracilResidue):
    ToReturn = {}
    UracilConsecutiveAtoms = {"N1": ["C6", "C2", "C1'"],
                              "C2": ["O2", "N3", "N1"],
                              "O2": ["C2"],
                              "N3": ["H3", "C4", "C2"],
                              "H3": ["N3"],
                              "C4": ["O4", "C5", "N3"],
                              "O4": ["C4"],
                              "C5": ["H5", "C6", "C4"],
                              "H5": ["C5"],
                              "C6": ["H6", "N1", "C5"],
                              "H6": ["C6"],
                               "P": ["OP1", "O5'", "OP2"],
                               "OP1": ["P"],
                               "OP2": ["P"],
                               "O5'": ["P", "C5'"],
                               "C5'": ["H5'", "H5''", "C4'", "O5'"],
                               "C4'": ["H4'", "C3'", "O4'", "C5'"],
                               "C3'": ["O3'", "H3'", "C2'", "C4'"],
                               "O3'": ["C3'"],
                               "O4'": ["C4'", "C1'"],
                               "C2'": ["H2'", "C1'", "O2'", "C3'"],
                               "O2'": ["HO2'", "C2'"],
                               "H1'": ["C1'"],
                               "H2'": ["C2'"],
                               "H3'": ["C3'"],
                               "H4'": ["C4'"],
                               "H5'": ["C5'"],
                               "H5''": ["C5'"],
                               "HO2'": ["O2'"],
                               "C1'": ["O4'", "H1'", "C2'", "N1"]}

    # Find out out HO2' is
    # Add hydrogens degree 1
    for atom in UracilResidue:
        AtomNameClean = atom.get_fullname().replace(" ", "").replace("''", "'")
        if AtomNameClean in UracilConsecutiveAtoms:
            AtomsBondedToCurrentAtom = UracilConsecutiveAtoms[AtomNameClean]
            VerticesToAdd = []
            for i in range(len(AtomsBondedToCurrentAtom)):
                try:
                    VerticesToAdd.append(UracilResidue[AtomsBondedToCurrentAtom[i]])
                except:
                    pass
            ToReturn[atom] = VerticesToAdd
    return ToReturn

def BuildMoleculeGraph(pdbReader):
    ResiNum = 1
    ResidueList = []
    for model in structure:
        for chain in model:  # all da chains
            for residue in chain:  # all da resis in da chain
                ResiName = "{}{}".format(ResiNum, residue.get_resname())
                if residue.get_resname() == "G":
                    ResidueList.append(ResiGraph(ResiName, ConstructGuanineConsecutive(residue)))

                if residue.get_resname() == "A":
                    ResidueList.append(ResiGraph(ResiName, ConstructAdenineConsecutive(residue)))

                if residue.get_resname() == "C":
                    ResidueList.append(ResiGraph(ResiName, ConstructCytosineConsecutive(residue)))

                if residue.get_resname() == "U":
                    ResidueList.append(ResiGraph(ResiName, ConstructUracilConsecutive(residue)))
                ResiNum += 1
        break

    return ResidueList

def FormatGraphsForC(CurrentResidue):
    Keys = list(CurrentResidue.AtomLocationsInAdjMatrix.keys())
    for key in Keys:
        Neighbors = list(CurrentResidue.AtomLocationsInAdjMatrix[key])
        for neighbor in Neighbors:
            CurrentRootCoords = list(key.get_vector().get_array())
            CurrentNeighborCoords = list(neighbor.get_vector().get_array())
            print("{},{},{},{},{},{},{},{},{}".format(CurrentResidue.Name, key.get_fullname().strip(), CurrentRootCoords[0], CurrentRootCoords[1], CurrentRootCoords[2], neighbor.get_fullname().strip(), CurrentNeighborCoords[0], CurrentNeighborCoords[1], CurrentNeighborCoords[2]))


if __name__ == '__main__':
    structure = PDBParser().get_structure("2n7x", r"C:\Users\Trist\Downloads\2n7x.pdb")
    ListOfResidueGraphs = BuildMoleculeGraph(structure)
    
    for Resi in ListOfResidueGraphs:
        FormatGraphsForC(Resi)
