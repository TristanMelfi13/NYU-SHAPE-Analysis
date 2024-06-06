using Microsoft.VisualBasic.FileIO;
using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Diagnostics.Tracing;
using System.Linq;
using System.Net.Http.Headers;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;
using static System.Net.Mime.MediaTypeNames;

namespace NYU_SHAPE_Analysis
{
    internal class NYUMain
    {
        static List<DictionarySet<string, double>> FinalDictionarySets = new();

        static Dictionary<string, LinkedList<Residue>> ListOfMolecules = new();

        static Dictionary<string, Dictionary<string, ExperimentalData>> ListOfExperimentalData = new();

        static Dictionary<string, List<Tuple<string, string>>> NonCanonicalPairs = new();

        static Dictionary<string, List<Tuple<string, string>>> WCPairs = new();

        static Dictionary<string, List<string>> HairpinNucleotides = new();

        static void AddToMolecules(string Name, LinkedList<Residue> residues)
        {
            lock (ListOfMolecules)
            {
                ListOfMolecules.Add(Name, residues);
            }
        }


        static Dictionary<string, int> ResiNumToResi = new Dictionary<string, int>();
        static Dictionary<TKey, TValue> SortDictionaryByValue<TKey, TValue>(Dictionary<TKey, TValue> dictionary)
        {
            return dictionary.OrderBy(pair => pair.Value).ToDictionary(pair => pair.Key, pair => pair.Value);
        }
        static int GetJusResiNumber(string str)
        {
            LinkedList<string> ToReturn = new LinkedList<string>();
            foreach (char c in str)
            {
                if (char.IsLetter(c))
                {
                    return int.Parse(string.Join("", ToReturn));
                }
                else
                {
                    ToReturn.AddLast(c.ToString());
                }
            }
            throw new Exception("This should never have failed what did you input?");
        }
        static void GetFullMoleculeGraph(string Molecule, string MoleculeName)
        {
            string pythonExePath = "python.exe"; // Set the correct path
            string scriptPath = "BioPythonPDBGetter.py"; // Set the path to your Python script

            ProcessStartInfo startInfo = new ProcessStartInfo
            {
                FileName = pythonExePath,
                Arguments = $"{scriptPath} {Molecule}",
                UseShellExecute = false,
                RedirectStandardOutput = true
            };

            using (Process process = new Process { StartInfo = startInfo })
            {
                process.Start();
                string output = process.StandardOutput.ReadToEnd();
                var ResiGraph = BuildGraphs(output);
                AddToMolecules(MoleculeName, ResiGraph);

            }
        }
        static void GetReactivityData(string MoleculeName)
        {
            string pythonExePath = "python.exe"; // Set the correct path
            string scriptPath = "ReactivityFileGetter.py"; // Set the path to your Python script


            MoleculeName = MoleculeName.Substring(0, MoleculeName.IndexOf("."));

            string ReactivityFilePath = @"DataSets\" + MoleculeName + "_Reactivity_Data.csv";


            ProcessStartInfo startInfo = new ProcessStartInfo
            {
                FileName = pythonExePath,
                Arguments = $"{scriptPath} {ReactivityFilePath}",
                UseShellExecute = false,
                RedirectStandardOutput = true
            };

            using (Process process = new Process { StartInfo = startInfo })
            {
                process.Start();
                string output = process.StandardOutput.ReadToEnd();

                var ExpData = new LinkedList<ExperimentalData>(); // A list of ExperimentalData?
                string[] data = output.Split("*");

                foreach (string d in data)
                {
                    if (d.Length != 0)
                    {
                        ExperimentalData NewData = CreateDataClass(d);
                        NewData.CreateAvgDatas();
                        ExpData.AddLast(NewData);
                        if (ListOfExperimentalData.ContainsKey(MoleculeName))
                        {
                            ListOfExperimentalData[MoleculeName].Add(NewData.GetReagentName(), NewData);

                        } else
                        {
                            ListOfExperimentalData.Add(MoleculeName, new Dictionary<string, ExperimentalData>());
                            ListOfExperimentalData[MoleculeName].Add(NewData.GetReagentName(), NewData);
                        }
                    }
                }
            }
        }
        static LinkedList<Residue> BuildGraphs(string OutputFromPython)
        {
            string[] Lines = OutputFromPython.Split(new char[] { '\n', });
            LinkedList<Residue> ListOfResidues = new();
            string CurrentResidueName = "";
            Residue residue = new Residue(CurrentResidueName);
            foreach (string line in Lines)
            {
                string[] Info = line.Split(new char[] { ',' });
                string ToCompare = Info[0];

                if (Info[0].Length == 0)
                {
                    break;
                }
                if (!CurrentResidueName.Equals(ToCompare))
                {
                    CurrentResidueName = ToCompare;
                    // Add the Residue We have been building
                    ListOfResidues.AddLast(residue);
                    // Create a new Residue
                    residue = new Residue(ToCompare);
                }
                residue.AddAtom(Info);
            }
            ListOfResidues.AddLast(residue);


            ListOfResidues.RemoveFirst();
            return ListOfResidues;
        }
        static void DoMath(Residue CurrentResidue, LinkedList<Residue> ListOfAllResidues, string Path)
        {
            Func<double[], double[], double> Distance = (v1, v2) => { return Math.Sqrt(Math.Pow(v1[0] - v2[0], 2) + Math.Pow(v1[1] - v2[1], 2) + Math.Pow(v1[2] - v2[2], 2)); };
            string CurrentResidueName = CurrentResidue.GetName();
            LinkedList<Atom> CurrentResidueAtoms = CurrentResidue.GetAtomList();
            Path += "\\" + CurrentResidue.GetName();
            System.IO.Directory.CreateDirectory(Path);
            // First Things First Calculate allllll Da Intranucleotide lengths
            foreach (Atom atom1 in CurrentResidueAtoms)
            {
                string CurrentFile = Path + "\\" + atom1.GetAtomName() + ".txt";
                LinkedList<string> Data = new LinkedList<string>();
                foreach (Residue residue in ListOfAllResidues)
                {
                    foreach (Atom atom2 in residue.GetAtomList())
                    {
                        double dist = Distance(atom1.GetAtomCoords(), atom2.GetAtomCoords());
                        string MeowMeow = CurrentResidue.GetName() + "\t" + atom1.GetAtomName() + "\t" + residue.GetName() + "\t" + atom2.GetAtomName() + "\t" + dist + "\n";
                        Data.AddLast(MeowMeow);
                    }
                }
                System.IO.File.WriteAllText(CurrentFile, string.Join("\n", Data));
            }
        }

        static void ParralelDoMath(LinkedList<Residue> Molecule, string PathName)
        {

            foreach (var resi in Molecule)
            {
                DoMath(resi, Molecule, PathName);
            }
        }

        static string[] SubSplit(string ThingToSplit, char Delimeter)
        {
            return ThingToSplit.Split(new char[] { Delimeter });
        }
        static ExperimentalData CreateDataClass(string DataForReaction)
        {
            ExperimentalData ToReturn = new ExperimentalData(DataForReaction.Substring(0, DataForReaction.IndexOf("\n")));
            string[] ActualValues = DataForReaction.Split("%");
            for (int i = 1; i < ActualValues.Length; i++)
            {
                string str = string.Join(" ", SubSplit(ActualValues[i], '%'));
                string[] AlmostCleaned = SubSplit(str, '\n');

                // Concentration is AlmostCleaned[0]

                string Concentration = AlmostCleaned[0];


                LinkedList<Tuple<string, double>> Results = new LinkedList<Tuple<string, double>>();
                for (int j = 1; j < AlmostCleaned.Length; j++)
                {
                    string[] FinalSplit = SubSplit(AlmostCleaned[j], ',');
                    try
                    {
                        var NucleotideValuePair = new Tuple<string, double>(FinalSplit[0], double.Parse(FinalSplit[1]));
                        Results.AddLast(NucleotideValuePair);
                    }
                    catch (Exception)
                    {

                    }
                }

                if (Concentration.Contains("nan"))
                {
                    return ToReturn;
                }
                ToReturn.AddToResults(Concentration.Substring(Concentration.IndexOf("_") + 1), Results);
            }

            return ToReturn;
        }
        static Dictionary<string, double> GetSHAPEData(LinkedList<ExperimentalData> Guys)
        {
            Dictionary<string, double> ToReturn = new Dictionary<string, double>();
            Func<double[], int, double> Avg = (v1, v2) =>
            {
                double RunningAvg = 0;
                foreach (double d in v1)
                {
                    RunningAvg += d;
                }
                return RunningAvg / v2;
            };
            foreach (var kvp in Guys.First().GetAvgVals())
            {
                // For every nucleotide in our molecule

                string CurrentNucleotide = kvp.Key;

                double RunninSum = 0;
                LinkedList<double> ForPrinting = new LinkedList<double>();
                foreach (ExperimentalData dude in Guys)
                {

                    try
                    {
                        RunninSum += dude.GetAvgVals()[CurrentNucleotide];
                        ForPrinting.AddLast(dude.GetAvgVals()[CurrentNucleotide]);
                    }
                    catch (Exception)
                    {

                    }
                }
                double Avgerage = RunninSum / Guys.Count();
                ToReturn.Add(CurrentNucleotide, Avgerage);
                // Console.WriteLine(CurrentNucleotide + " has an avg reactivity value ---> " + Avgerage);
            }


            return ToReturn;
        }
        static void RemoveKeysFromDictionary(Dictionary<string, double> dictionary, params string[] keysToRemove)
        {
            foreach (var key in keysToRemove)
            {
                dictionary.Remove(key);
            }
        }
        static string GetNucleotideType(string input)
        {
            if (input.Contains("A"))
            {
                return "U";
            }
            else if (input.Contains("U"))
            {
                return "A";
            }
            else if (input.Contains("G"))
            {
                return "C";
            }
            else
            {
                return "G";
            }
        }
        static string GetSelfNucleotideType(string input)
        {
            if (input.Contains("A"))
            {
                return "A";
            }
            else if (input.Contains("U"))
            {
                return "U";
            }
            else if (input.Contains("G"))
            {
                return "G";
            }
            else
            {
                return "C";
            }
        }
        static string[] FindMostLikelyPair(string CurrentNucleotide, List<string> PotentialPairs, int SequenceCount)
        {
            if (PotentialPairs.Count() == 0)
            {
                return new string[] { CurrentNucleotide, "meow" };
            }
            int CurrentResiNumber = GetJusResiNumber(CurrentNucleotide); // Try to match him best as possible
            var ForMatching = new Dictionary<string, int>();
            foreach (var resi in PotentialPairs)
            {
                int RelativeSequencePosition = SequenceCount - GetJusResiNumber(resi);
                int Difference = Math.Abs(CurrentResiNumber - RelativeSequencePosition);
                if (!ForMatching.ContainsKey(resi))
                {
                    ForMatching.Add(resi, Difference);
                }
            }
            ForMatching = SortDictionaryByValue(ForMatching);
            foreach (var kvp in ForMatching)
            {
                if (kvp.Value != 0)
                {
                    return new string[] { CurrentNucleotide, "meow" };
                }
                else
                {
                    return new string[] { CurrentNucleotide, kvp.Key.ToString() };

                }
            }
            throw new Exception("This was not supposed to happen... ");
        }
        static Dictionary<string, double> FilterDictionaryByKey(List<string> ToMatch, Dictionary<string, double> DictionaryToFilterOn)
        {
            var ToReturn = new Dictionary<string, double>();
            foreach (string str in ToMatch)
            {
                if (DictionaryToFilterOn.ContainsKey(str))
                {
                    ToReturn.Add(str, DictionaryToFilterOn[str]);
                }
            }
            return ToReturn;
        }
        static Queue<string> RemoveSettledNucleotides(Queue<string> queue, string[] PairsToRemove)
        {
            var ToRemove = new Queue<string>();

            foreach (var nucleotide in queue)
            {
                if (!nucleotide.Equals(PairsToRemove[0]) && !nucleotide.Equals(PairsToRemove[1]))
                {
                    ToRemove.Enqueue(nucleotide);
                }
            }
            return ToRemove;
        }

        static List<Tuple<string, string>> GetBonded(List<string> NucleotidesToMatch)
        {
            var MakeSureNoDoubles = new HashSet<string>();
            var ToReturn = new List<Tuple<string, string>>();

            foreach (var Nucleotide in NucleotidesToMatch)
            {
                foreach (var Nucleotide2 in NucleotidesToMatch)
                {
                    var NucleotideNumber = GetJusResiNumber(Nucleotide);
                    var Nucleotide2Number = GetJusResiNumber(Nucleotide2);


                    if (!MakeSureNoDoubles.Contains(Nucleotide + Nucleotide2))
                    {
                        if (NucleotideNumber == (Nucleotide2Number + 1) || NucleotideNumber == (Nucleotide2Number - 1))
                        {
                            MakeSureNoDoubles.Add(Nucleotide + Nucleotide2);
                            MakeSureNoDoubles.Add(Nucleotide2 + Nucleotide);
                            ToReturn.Add(new Tuple<string, string>(Nucleotide, Nucleotide2));
                        }
                    }
                }
            }
            return ToReturn;
        }

        static void ComputeInteractions(Dictionary<string, string> Directories, List<Tuple<string, string>> ResiduesComputingOn, DictionarySet<string, double> DSet, int SequenceCount, ConsolidatedResults Results)
        {
            foreach (var Pair in ResiduesComputingOn)
            {
                var SetOfAdded = new HashSet<string>();

                var Item1 = Pair.Item1;
                var Item2 = Pair.Item2;


                if (GetJusResiNumber(Pair.Item1) == 1 || GetJusResiNumber(Pair.Item1) == SequenceCount)
                {
                    var temp = Item2;
                    Item2 = Item1;
                    Item1 = temp;
                }

                var DirectoryForInteractions = Directories[Item1];
                var AllFiles = Directory.GetFiles(DirectoryForInteractions);

                foreach (var file in AllFiles)
                {
                    var FileContent = File.ReadAllLines(file);
                    foreach (var line in FileContent)
                    {
                        if (line.Length != 0)
                        {
                            var split = line.Split("\t");

                            if (split[2].Equals(Item2))
                            {

                                if (!SetOfAdded.Contains(split[1] + split[3]) && !SetOfAdded.Contains(split[3] + split[1]))
                                {
                                    SetOfAdded.Add(split[1] + split[3]);
                                    DSet.AddToDictionary(split[0] + "\t" + split[2], split[1] + split[3], double.Parse(split[4]));

                                }
                            }
                        }
                    }
                }
                SetOfAdded.Clear();
            }

            Results.AddToResults(DSet);

        }

        static Dictionary<string, string> NucleotidesToPath(string[] directories)
        {
            var ToReturn = new Dictionary<string, string>();


            foreach (var d in directories)
            {
                ToReturn.Add(d.Substring(d.LastIndexOf("\\") + 1), d);
            }
            return ToReturn;
        }

        static string FindWCNucleotide(string CurrentNucleotide, HashSet<string> ListToSearchThrough, int SequenceCount)
        {
            int NucleotideNumber = GetJusResiNumber(CurrentNucleotide);
            int OtherSideNucleotide = Math.Abs(NucleotideNumber - (SequenceCount));
            if (CurrentNucleotide.Contains("G") && ListToSearchThrough.Contains(OtherSideNucleotide + "C"))
            {
                return OtherSideNucleotide + "C";
            }

            if (CurrentNucleotide.Contains("C") && ListToSearchThrough.Contains(OtherSideNucleotide + "G"))
            {
                return OtherSideNucleotide + "G";
            }

            if (CurrentNucleotide.Contains("A") && ListToSearchThrough.Contains(OtherSideNucleotide + "U"))
            {
                return OtherSideNucleotide + "U";
            }

            if (CurrentNucleotide.Contains("U") && ListToSearchThrough.Contains(OtherSideNucleotide + "A"))
            {
                return OtherSideNucleotide + "A";
            }

            return "null";

        }

        static List<Tuple<string, string>> FindWatsonCrickPairs(List<string> UnreactiveNucleotides, List<string> ReactiveNucleotides)
        {
            var ToReturn = new List<Tuple<string, string>>();
            // Check for WC in Reactive -=-=-> Unreactive
            var UnreactiveHashSet = UnreactiveNucleotides.ToHashSet();
            var CopyToWorkOn = UnreactiveNucleotides;
            int SeqCount = UnreactiveNucleotides.Count() + ReactiveNucleotides.Count() + 1;
            
            
            foreach (var Nucleotide in ReactiveNucleotides)
            {
                var CurrentResiNumber = GetJusResiNumber(Nucleotide);
                var OtherWCNucleotide = FindWCNucleotide(Nucleotide, UnreactiveHashSet, SeqCount);
                if (!OtherWCNucleotide.Equals("null"))
                {
                    ToReturn.Add(new Tuple<string, string>(Nucleotide, OtherWCNucleotide));
                    UnreactiveHashSet.Remove(OtherWCNucleotide);
                }
            }
            foreach (var Nucleotide in UnreactiveNucleotides)
            {
                var WCNucleotide = FindWCNucleotide(Nucleotide, UnreactiveHashSet, SeqCount);
                var Perm1 = new Tuple<string, string>(Nucleotide, WCNucleotide);
                var Perm2 = new Tuple<string, string>(WCNucleotide, Nucleotide);
                if (!WCNucleotide.Equals("null") && (!ToReturn.Contains(Perm1) && !ToReturn.Contains(Perm2)))
                {
                    ToReturn.Add(Perm1);
                }
            }
            return ToReturn;
        }

        static Tuple<string, string> FindTerminalNucleotides(Tuple<List<string>, List<string>> UnreactiveAndReactiveNucleotides, int SeqCount)
        {
            var Unreactive = UnreactiveAndReactiveNucleotides.Item1;
            var Reactive = UnreactiveAndReactiveNucleotides.Item2;

            var FirstNucleotide = "";

            if (Unreactive.First().Contains("1"))
            {
                FirstNucleotide = Unreactive.First();
            } else if (Reactive.First().Contains("1"))
            {
                FirstNucleotide = Reactive.First();
            } else
            {
                throw new Exception("Something went wrong in FindTerminalNucleotides");
            }

            var LastNucleotide = "";
            if (Unreactive.Last().Contains(SeqCount + "" + GetNucleotideType(FirstNucleotide)))
            {
                LastNucleotide = Unreactive.Last();
            } else if (Reactive.Last().Contains(SeqCount + "" + GetNucleotideType(FirstNucleotide)))
            {
                LastNucleotide = Reactive.Last();
            }
            return new Tuple<string, string>(FirstNucleotide, LastNucleotide);
        }

        static Tuple<string, string> FindMatchingTuple(string MissingHisFriend, List<Tuple<string, string>> ToSearchThrough)
        {
            foreach(var pair in ToSearchThrough)
            {
                if (pair.Item1.Equals(MissingHisFriend) || pair.Item2.Equals(MissingHisFriend))
                {
                    return pair;
                }
            }

            throw new Exception("No mathcing tuple");


        }


        static Tuple<string, string> GetLastWCBeforeHairpin(List<Tuple<string, string>> WCPairs, int SeqCount)
        {
            var Middle = SeqCount / 2;
            var ForLookupLaterPos = new Dictionary<string, int>();
            var ForLookupLaterNeg = new Dictionary<string, int>();
            foreach (var Pair in WCPairs)
            {
                var Diff1 = GetJusResiNumber(Pair.Item1) - Middle;
                var Diff2 = GetJusResiNumber(Pair.Item2) - Middle;
                if (Diff1 > 0)
                {
                    ForLookupLaterPos.Add(Pair.Item1, Diff1);
                } else
                {
                    ForLookupLaterNeg.Add(Pair.Item1, Diff1);
                }

                if (Diff2 > 0)
                {
                    ForLookupLaterPos.Add(Pair.Item2, Diff2);
                }
                else
                {
                    ForLookupLaterNeg.Add(Pair.Item2, Diff2);
                }
            }
            SortDictionaryByValue(ForLookupLaterNeg);
            var LastWCNucleotide1 = "";
            foreach(var kvp in ForLookupLaterNeg)
            {
                LastWCNucleotide1 = kvp.Key;
                break;
            }
            return FindMatchingTuple(LastWCNucleotide1, WCPairs);
        }


        static List<string> FindHairpin(List<Tuple<string, string>> WCPairs, Tuple<string, string> TerminalNucleotides, int SeqCount, Tuple<List<string>, List<string>> UnreactiveAndReactiveNucleotides)
        {
            var CopyToWorkWith = UnreactiveAndReactiveNucleotides.Item2;
            var LastWC = GetLastWCBeforeHairpin(WCPairs, SeqCount);
            var WCCopy = WCPairs;
            WCCopy.Add(TerminalNucleotides);
            foreach(var Pair in WCCopy)
            {
                if (CopyToWorkWith.Contains(Pair.Item1))
                {
                    CopyToWorkWith.Remove(Pair.Item1);
                }
                if (CopyToWorkWith.Contains(Pair.Item2))
                {
                    CopyToWorkWith.Remove(Pair.Item2);
                }
            }
            CopyToWorkWith.Add(LastWC.Item1);
            CopyToWorkWith.Add(LastWC.Item2);
            CopyToWorkWith = SortNucleotides(CopyToWorkWith);
            var UnreactiveCopyArray = CopyToWorkWith.ToArray();
            for (int i = 0; i < UnreactiveCopyArray.Length - 1; i++)
            {
                var Diff = GetJusResiNumber(UnreactiveCopyArray[i + 1]) - GetJusResiNumber(UnreactiveCopyArray[i]);

                if (Diff >= 4)
                {
                    CopyToWorkWith.Remove(UnreactiveCopyArray[i]);
                }
                if (i == UnreactiveCopyArray.Length - 2 && Diff >= 4)
                {
                    CopyToWorkWith.Remove(UnreactiveCopyArray[i + 1]);
                }
            }
            FillInHairpinGaps(CopyToWorkWith, UnreactiveAndReactiveNucleotides.Item1);
            if (!CopyToWorkWith.Contains(LastWC.Item1))
            {
                CopyToWorkWith.Add(LastWC.Item1);
            }

            if (!CopyToWorkWith.Contains(LastWC.Item2))
            {
                CopyToWorkWith.Add(LastWC.Item2);
            }

            return SortNucleotides(CopyToWorkWith);



        }

        static void FillInHairpinGaps(List<string> CurrentHairpinNucleotides, List<string> RemainingPotentialNucleotides)
        {
            var CurrentHairpinNucleotidesArray = CurrentHairpinNucleotides.ToArray();

            var RemainingPotentialNucleotidesHashSet = RemainingPotentialNucleotides.ToHashSet();

            for (int i = 0; i < CurrentHairpinNucleotidesArray.Length - 1; i++)
            {
                var Diff = GetJusResiNumber(CurrentHairpinNucleotidesArray[i + 1]) - GetJusResiNumber(CurrentHairpinNucleotidesArray[i]);

                if (Diff != 1)
                {

                    var CurrentResiNumber = GetJusResiNumber(CurrentHairpinNucleotidesArray[i]) + 1;

                    if (RemainingPotentialNucleotidesHashSet.Contains(CurrentResiNumber + "A"))
                    {
                        CurrentHairpinNucleotides.Insert(i, CurrentResiNumber + "A");

                    } else if (RemainingPotentialNucleotidesHashSet.Contains(CurrentResiNumber + "G"))
                    {
                        CurrentHairpinNucleotides.Insert(i, CurrentResiNumber + "G");

                    }
                    else if (RemainingPotentialNucleotidesHashSet.Contains(CurrentResiNumber + "C"))
                    {
                        CurrentHairpinNucleotides.Insert(i, CurrentResiNumber + "C");

                    }
                    else if (RemainingPotentialNucleotidesHashSet.Contains(CurrentResiNumber + "U"))
                    {
                        CurrentHairpinNucleotides.Insert(i, CurrentResiNumber + "U");
                    }
                }
            }
        }

        


        static void Compute(Tuple<List<string>, List<string>> UnreactiveAndReactiveNucleotides, string MoleculeName, string RName) // Each Reagent gets their own results...
        {
            int SeqCount = UnreactiveAndReactiveNucleotides.Item1.Count() + UnreactiveAndReactiveNucleotides.Item2.Count() ;
            var WCPairs = FindWatsonCrickPairs(UnreactiveAndReactiveNucleotides.Item1, UnreactiveAndReactiveNucleotides.Item2);
            var TerminalNucleotides = FindTerminalNucleotides(UnreactiveAndReactiveNucleotides, SeqCount);
            var HairpinNucleotides = FindHairpin(WCPairs, TerminalNucleotides, SeqCount, UnreactiveAndReactiveNucleotides);

            Console.WriteLine("Unreactive Nucleotides: " + string.Join(", ", UnreactiveAndReactiveNucleotides.Item1));
            Console.WriteLine("Reactive Nucleotides: " + string.Join(", ", UnreactiveAndReactiveNucleotides.Item2));
            Console.WriteLine("WCPairs: " + string.Join(", ", WCPairs));
            Console.WriteLine("Hairpin: " + string.Join(", ", HairpinNucleotides));



        }

        static void ExecuteThreads(List<Thread> ToStart)
        {
            foreach (var t in ToStart)
            {
                t.Start();
            }

            foreach (var t in ToStart)
            {
                t.Join();
            }

            ToStart.Clear();
        }
        static double Average(List<double> ToAverage)
        {
            double RunningSum = 0;
            foreach (var d in ToAverage)
            {
                RunningSum += d;
            }
            return RunningSum / ToAverage.Count();
        }
        
        static void CreateDirectory(string directory)
        {
            if (!Directory.Exists(directory))
            {
                Directory.CreateDirectory(directory);
            }
        }


        static void SaveResults(string OutputDirectory, Dictionary<string, ConsolidatedResults> Results)
        {
            // After we combine all the results this is the very last thing to do

            CreateDirectory(OutputDirectory);

            foreach (var result in Results)
            {
                var DirectoryToCreate = OutputDirectory + "\\" + result.Key;
                CreateDirectory(DirectoryToCreate);
                var CurrentResults = result.Value;
                var Keys = CurrentResults.GetDictionarySet().GetSuperKeys();
                foreach (var k in Keys)
                {
                    var SubSet = CurrentResults.GetDictionarySet().GetSubSet(k);
     
                    string FilePath = DirectoryToCreate + "\\" + k + ".txt";
                    var ToAdd = new List<string>();
                    foreach (var kvp in SubSet)
                    {
                        ToAdd.Add(k + "\t" + kvp.Key + "\t" + Average(kvp.Value));
                    }
                    File.WriteAllLines(FilePath, ToAdd);
                    ToAdd.Clear();
                }
            }
        }

        static Dictionary<string, string> PopulateMoleculeList()
        {
            var ToReturn = new Dictionary<string, string>();
            var Files = Directory.GetFiles("PDBFiles", "*.pdb");

            foreach (string f in Files)
            {
                string MoleculeName = f.Substring(f.IndexOf("\\") + 1);
                ToReturn.Add(MoleculeName, f);
            }
            return ToReturn;
        }

        static double StdDev(List<double> input, double avg)
        {
            double ToReturn = 0.0;
            foreach (double d in input)
            {
                ToReturn = Math.Pow(d - avg, 2);
            }
            ToReturn /= input.Count();
            return Math.Sqrt(ToReturn);


        }

        static void DoubleCheck(Tuple<List<string>, List<string>> Nucleotides, string MoleculeName)
        {
            Console.WriteLine(MoleculeName);
            int SeqCount = Nucleotides.Item1.Count() + Nucleotides.Item2.Count();
            var HalfSeqCount = Math.Ceiling(SeqCount / 2.0);

            /*Console.WriteLine(string.Join(", ", Nucleotides.Item1));
            Console.WriteLine(string.Join(", ", Nucleotides.Item2));*/
        }


        static List<string> SortNucleotides(List<string> ToSort)
        {
            return ToSort.OrderBy(item => int.Parse(new String(item.Where(Char.IsDigit).ToArray()))).ThenBy(item => item.Where(Char.IsLetter).ToArray()).ToList();
        }
        static double GetJumpLimit(double[] ToComputeOn)
        {
            List<double> Differences = new List<double>();

            for (int i = 0; i < ToComputeOn.Length - 1; i++)
            {
                Differences.Add(ToComputeOn[i + 1] - ToComputeOn[i]);
            }

            var avg = Average(Differences);

            return avg;
        }      
        static Tuple<List<string>, List<string>> DetermineLowReactivityNucleotides(Dictionary<string, double> NucleotideReactivityPairs, string MoleculeName)
        {
            var ReactivityData = SortDictionaryByValue(NucleotideReactivityPairs).Values.ToArray();
            var ReactivityDataKeys = SortDictionaryByValue(NucleotideReactivityPairs).Keys.ToArray();
            HashSet<double> ForAnalysisLater = new HashSet<double>();
            var Avg = GetJumpLimit(ReactivityData.ToArray());
            var LowReactivityNucleotides = new List<string>();
            var HighReactivityNucleotides = new List<string>();
            for (int i = 0; i < ReactivityData.Length - 1; i++)
            {
                LowReactivityNucleotides.Add(ReactivityDataKeys.ElementAt(i));
                var Difference = ReactivityData[i + 1] - ReactivityData[i];
                if (Difference >= Math.Round(Avg, 1))
                {
                    // HighReactivityNucleotides.Add(ReactivityDataKeys.ElementAt(i));
                    for (int j = i + 1; j < ReactivityDataKeys.Length; j++)
                    {
                        HighReactivityNucleotides.Add(ReactivityDataKeys.ElementAt(j));
                    }
                    break;
                }
            }
            // HighReactivityNucleotides.RemoveAt(HighReactivityNucleotides.Count() - 1);
            return new Tuple<List<string>, List<string>>(SortNucleotides(LowReactivityNucleotides), SortNucleotides(HighReactivityNucleotides));
        }








        static Dictionary<string, string> NucleotideToPath(List<string> Directories)
        {
            var ToReturn = new Dictionary<string, string>();


            foreach (var s in Directories)
            {
                ToReturn.Add(s.Substring(s.LastIndexOf("\\") + 1), s);
            }

            return ToReturn;


        }


        static void Main()
        {
            
            /*Populates ListOfMolecules with kvp's name to a LinkedList of Nucleotides*/

            foreach(var kvp in PopulateMoleculeList())
            {
                GetFullMoleculeGraph(kvp.Value, kvp.Key); 
            }


            /*List<Thread> Workers = new List<Thread>();

            foreach(var Molecule in ListOfMolecules)
            {
                Workers.Add(new Thread(() => ParralelDoMath(Molecule.Value, "DistanceOutputs\\" + Molecule.Key.Substring(0, Molecule.Key.IndexOf(".")))));
            }
            ExecuteThreads(Workers);*/


            foreach(var kvp in ListOfMolecules)
            {
                GetReactivityData(kvp.Key);
            }

            foreach(var Molecule in ListOfExperimentalData) // For every molecule we have a file for
            {
                Console.WriteLine(Molecule.Key);
                foreach (var Reagent in Molecule.Value) // For every entry in the molecule file
                {
                    var NucleotideReactivityPairs = Reagent.Value.GetAvgVals();

                    // Console.WriteLine(string.Join("\n", NucleotideReactivityPairs));

                    var LowAndHighReactivityNucleotides = DetermineLowReactivityNucleotides(NucleotideReactivityPairs, Molecule.Key);
                    var LowReactivity = LowAndHighReactivityNucleotides.Item1;
                    var HighReactivity = LowAndHighReactivityNucleotides.Item2;
                    Compute(LowAndHighReactivityNucleotides, Molecule.Key, Reagent.Key.Replace("\r", ""));
                    break;
                }
                break;
            }
            
        }
    }
}
