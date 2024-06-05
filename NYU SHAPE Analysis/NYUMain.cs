using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class NYUMain
    {

        static string GetFullMoleculeGraph()
        {
            string pythonExePath = "python.exe"; // Set the correct path
            string scriptPath = "ConsecutiveAtomPairs.py"; // Set the path to your Python script

            ProcessStartInfo startInfo = new ProcessStartInfo
            {
                FileName = pythonExePath,
                Arguments = scriptPath,
                UseShellExecute = false,
                RedirectStandardOutput = true
            };

            using (Process process = new Process { StartInfo = startInfo })
            {
                process.Start();
                string output = process.StandardOutput.ReadToEnd();
                return output;
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
            ListOfResidues.RemoveFirst();
            return ListOfResidues;
        }

        static void DoMath(Residue CurrentResidue, LinkedList<Residue> ListOfAllResidues)
        {
            Func<double[], double[], double> Distance = (v1, v2) => { return Math.Sqrt(Math.Pow(v1[0] - v2[0], 2) + Math.Pow(v1[1] - v2[1], 2) + Math.Pow(v1[2] - v2[2], 2)); };
            string CurrentResidueName = CurrentResidue.GetName();
            LinkedList<Atom> CurrentResidueAtoms = CurrentResidue.GetAtomList();
            string Path = @"Outputs\" + CurrentResidueName;
            System.IO.Directory.CreateDirectory(Path);

            // First Things First Calculate allllll Da Intranucleotide lengths
            foreach (Atom atom1 in CurrentResidueAtoms)
            {
                string CurrentFile = Path + "\\" + atom1.GetAtomName();
                foreach(Residue residue in ListOfAllResidues)
                {
                    foreach(Atom atom2 in residue.GetAtomList())
                    {
                        double dist = Distance(atom1.GetAtomCoords(), atom2.GetAtomCoords());
                    }
                }
            }








        }

        static void Main()
        {
            Stopwatch sw = Stopwatch.StartNew();
            string FullGraph = GetFullMoleculeGraph();
            string[] Lines = FullGraph.Split(new char[] { '\n', });
            LinkedList<Residue> ListOfResidues = BuildGraphs(FullGraph);
            LinkedList<Thread> Workers = new LinkedList<Thread>();
            System.IO.Directory.CreateDirectory("Outputs");



            foreach (Residue r in ListOfResidues)
            {
                Thread NewThread = new Thread(() => DoMath(r, ListOfResidues));
                NewThread.Name = r.GetName();
                Workers.AddLast(NewThread);
            }

            foreach(Thread thread in Workers)
            {
                thread.Start();
            }

            foreach(Thread thread in Workers)
            {
                thread.Join();
            }


            Console.WriteLine("Execution time: " + sw.ElapsedMilliseconds / 1000.0);
            // Make a thread for each residue

        }
    }
}
