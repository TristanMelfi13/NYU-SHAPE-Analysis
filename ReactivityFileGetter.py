import pandas as pd
reactivity_data = pd.read_csv("2N7X Reactivity Data.csv")

class ReactivityJunk:
    ReagentName = ""
    ReactionData = {}
    def __init__(self, reagentname):
        self.ReagentName = reagentname
    def AddData(self, Concentration, TupleList):
        self.ReactionData[Concentration] = TupleList


def CheckIfNumber(GuyToCheck):
    try:
        float(GuyToCheck)
        return True
    except:
        return False

def GetReagentNames(ColumnForSorting):
    Indices = [0]
    ToReturn = []
    count = 0
    for i in range(len(ColumnForSorting)):
        if not CheckIfNumber(ColumnForSorting[i]):
            if i == 0:
                pass
            if count % 2 == 0:
                ToReturn.append(ColumnForSorting[i])
            count += 1
    return ToReturn


def GetIndices(ColumnForSorting):
    Indices = [0]
    ToReturn = []
    for i in range(len(ColumnForSorting)):
        if not CheckIfNumber(ColumnForSorting[i]):
            if i == 0:
                pass
            elif i != Indices[-1] + 1:
                Indices.append(i)
    for i in range(len(Indices) - 1):
        ToReturn.append((Indices[i], Indices[i + 1] - 2))
    return ToReturn

def GetNucleotideNamesInSequence(DataFrame):
    last_column = DataFrame.iloc[:, len(DataFrame.columns) - 1].tolist()
    ToReturn = []
    for Guy in last_column:
        if Guy == "G" or Guy == "A" or Guy == "U" or Guy == "C":
            ToReturn.append(Guy)
    return ToReturn


def GetOnlyMeans(DataFrame):
    ColumnsToUse = []
    for i in range(0, len(DataFrame.columns), 2):
        SelectedRows = DataFrame.iloc[:, i].tolist()
        ColumnsToUse.append(SelectedRows)
    return ColumnsToUse

def ExtractReagentSpecificInfo(DataFrame, Nucleotides, ReagentName):
    ImportantInfo = DataFrame.iloc[0].tolist()
    ImportantInfo = ImportantInfo[::2]
    if ImportantInfo[-1] == "nucleotide":
        ImportantInfo.pop()
    DataFrame = DataFrame.iloc[1:] # Remove first row
    DataFrame = DataFrame.iloc[:, :-1] # Remove very last column which is usually missing nucleotides
    ReactivityValues = GetOnlyMeans(DataFrame)
    ToAddTo = ReactivityJunk(CurrentName)
    print("*{}".format(CurrentName))
    for i in range(len(ImportantInfo)):
        print("%{}".format(ImportantInfo[i]))
        TupleList = []  # Create a new TupleList for each iteration
        try:
            CurrentReactivities = ReactivityValues[i]
        except:
            break
        for j in range(len(CurrentReactivities)):
            print("{}{},{}".format(j + 1, Nucleotides[j], CurrentReactivities[j]))
            TupleList.append(("{}{}".format(j + 1, Nucleotides[j]), CurrentReactivities[j]))
        ToAddTo.AddData(ImportantInfo[i], TupleList)
    return ToAddTo


if __name__ == '__main__':
    reactivity_data = pd.read_csv("2N7X Reactivity Data.csv")
    cfs = list(reactivity_data.iloc[:, 0])
    cfs.insert(0, reactivity_data.columns[0])
    DataRanges = GetIndices(cfs)
    DataRanges.append((DataRanges[-1][1] + 2, reactivity_data.index[-1]))
    ReagentNames = GetReagentNames(cfs)
    NucleotideNames = GetNucleotideNamesInSequence(reactivity_data)
    Reactions = []
    for Tuple in DataRanges:
        CurrentName = ReagentNames.pop(0)
        Reactions.append(ExtractReagentSpecificInfo(reactivity_data.iloc[Tuple[0]:Tuple[1]], NucleotideNames, CurrentName))