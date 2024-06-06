using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class Residue
    {
        private string ResidueName = "";

        private LinkedList<Atom> AtomsInResidue = new LinkedList<Atom>();

        private Dictionary<string, Atom> KeyValuePairsAtoms = new Dictionary<string, Atom>();

        private Dictionary<string, double[]> BondWeights = new Dictionary<string, double[]>();


        public LinkedList<Atom> GetAtomList()
        {
            return AtomsInResidue;
        }
        
        public Residue(string residuename)
        {
            ResidueName = residuename;
        }

        public string GetName()
        {
            return ResidueName;
        }

        public Atom GetAtom(string AtomName) 
        { 
            return KeyValuePairsAtoms[AtomName];

        }

        public bool CheckPair(string Pair)
        {
            return BondWeights.ContainsKey(Pair);
        }


        public void AddAtom(string RootName, Atom InterResidueAtom, double[] DistanceBetweenThem)
        {
            var Atom1 = KeyValuePairsAtoms[RootName];
            Atom1.AddEdge(InterResidueAtom, DistanceBetweenThem);
            KeyValuePairsAtoms.Add(InterResidueAtom.GetAtomName() + "*", InterResidueAtom);
        }

        public void AddAtom(string RootName, string EdgeName, double[] DistanceBetweenThem)
        {



            Func<double[], double[]> Invert = (meow) => { return new double[] { -meow[0], -meow[1], -meow[2] }; };

            Atom? Atom1 = null;


            try
            {
                Atom1 = KeyValuePairsAtoms[RootName];
            }
            catch (Exception)
            {
                Atom NewAtom = new Atom(RootName, new double[] { 0, 0, 0 }, GetName());
                KeyValuePairsAtoms.Add(RootName, NewAtom);
                Atom1 = KeyValuePairsAtoms[RootName];
            }

            Atom? Atom2 = null;

            try
            {
                Atom2 = KeyValuePairsAtoms[EdgeName];
            } catch (Exception)
            {
                Atom NewAtom = new Atom(EdgeName, new double[] { 0, 0, 0}, GetName());
                KeyValuePairsAtoms.Add(EdgeName, NewAtom);
                Atom2 = KeyValuePairsAtoms[EdgeName];
            }


            if (!BondWeights.ContainsKey(Atom1.GetAtomName() + Atom2.GetAtomName()))
            {
                BondWeights.Add(Atom1.GetAtomName() + Atom2.GetAtomName(), DistanceBetweenThem);
            }






            Atom1.AddEdge(Atom2, DistanceBetweenThem);
        }


        public void AddAtom(string[] AtomicInformation)
        {
            string CurrentRootAtom = AtomicInformation[1];
            string CurrentNeighbor = AtomicInformation[5];

            Atom? atomRoot = null;
            Atom? atomNeighbor = null;
            foreach(Atom atom in AtomsInResidue)
            {
                if (atom.GetAtomName().Equals(CurrentRootAtom))
                {
                    atomRoot = atom;
                }

                if (atom.GetAtomName().Equals(CurrentNeighbor))
                {
                    atomNeighbor = atom;
                }

                if (atomRoot != null && atomNeighbor != null) 
                {
                    atomRoot.AddEdge(atomNeighbor);
                    return;
                }
            }


            if (atomRoot == null && atomNeighbor == null)
            {
                Atom NewRoot = new Atom(CurrentRootAtom, new double[] { double.Parse(AtomicInformation[2]), double.Parse(AtomicInformation[3]), double.Parse(AtomicInformation[4]) }, GetName());
                Atom NewNeighbor = new Atom(CurrentNeighbor, new double[] { double.Parse(AtomicInformation[6]), double.Parse(AtomicInformation[7]), double.Parse(AtomicInformation[8]) }, GetName());
                KeyValuePairsAtoms.Add(NewRoot.GetAtomName(), NewRoot);
                KeyValuePairsAtoms.Add(NewNeighbor.GetAtomName(), NewNeighbor);
                NewRoot.AddEdge(NewNeighbor);
                AtomsInResidue.AddLast(NewRoot);
                AtomsInResidue.AddLast(NewNeighbor);
                return;
            }
            else if (atomRoot != null && atomNeighbor == null)
            {
                Atom NewNeighbor = new Atom(CurrentNeighbor, new double[] { double.Parse(AtomicInformation[6]), double.Parse(AtomicInformation[7]), double.Parse(AtomicInformation[8]) }, GetName());
                KeyValuePairsAtoms.Add(NewNeighbor.GetAtomName(), NewNeighbor);

                atomRoot.AddEdge(NewNeighbor);
                AtomsInResidue.AddLast(NewNeighbor);
                return;

            }
            else if (atomRoot == null && atomNeighbor != null)
            {
                Atom NewRoot = new Atom(CurrentRootAtom, new double[] { double.Parse(AtomicInformation[2]), double.Parse(AtomicInformation[3]), double.Parse(AtomicInformation[4]) }, GetName());
                KeyValuePairsAtoms.Add(NewRoot.GetAtomName(), NewRoot);
                NewRoot.AddEdge(atomNeighbor);
                AtomsInResidue.AddLast(NewRoot);
                return;
            }
        }


        private Dictionary<string, LinkedList<double>> ComputeAverages(Atom[] CurrentBond, Dictionary<string, LinkedList<double>> ForComputing)
        {

            Func<double[], double[], double> Distance = (v1, v2) => { return Math.Sqrt(Math.Pow(v1[0] - v2[0], 2) + Math.Pow(v1[1] - v2[1], 2) + Math.Pow(v1[2] - v2[2], 2)); };

            if (ForComputing.ContainsKey(CurrentBond[0].GetAtomName() + CurrentBond[1].GetAtomName()) || ForComputing.ContainsKey(CurrentBond[1].GetAtomName() + CurrentBond[0].GetAtomName()))
            { 
                try
                {
                    ForComputing[CurrentBond[0].GetAtomName() + CurrentBond[1].GetAtomName()].AddLast(Distance(CurrentBond[0].GetRelativeAtomCoords(), CurrentBond[1].GetRelativeAtomCoords()));
                } catch (Exception) 
                {
                    ForComputing[CurrentBond[1].GetAtomName() + CurrentBond[0].GetAtomName()].AddLast(Distance(CurrentBond[1].GetRelativeAtomCoords(), CurrentBond[0].GetRelativeAtomCoords()));
                }
            } else
            {
                ForComputing.Add(CurrentBond[0].GetAtomName() + CurrentBond[1].GetAtomName(), new LinkedList<double> { });
                ForComputing[CurrentBond[0].GetAtomName() + CurrentBond[1].GetAtomName()].AddLast(Distance(CurrentBond[0].GetRelativeAtomCoords(), CurrentBond[1].GetRelativeAtomCoords()));
            }

            return ForComputing;

        }


        private double[] GetSpherical(double[] Cartesian)
        {
            double r = Math.Sqrt(Math.Pow(Cartesian[0], 2) + Math.Pow(Cartesian[1], 2) + Math.Pow(Cartesian[2], 2));
            double theta = Math.Atan2(Cartesian[1], Cartesian[0]);
            double phi = Math.Acos(Cartesian[2] / r);

            return new double[] { r, theta, phi };
        }




        public void PrintList(string ResiName)
        {


            Console.WriteLine(ResiName);

            if (ResidueName.Contains("C"))
            {
                Func<LinkedList<double>, double> Average = (ListOfDistances) =>
                {
                    double Avg = 0.0;
                    foreach (double d in ListOfDistances)
                    {
                        Avg += d;
                    }
                    return Avg / ListOfDistances.Count();
                };
                Dictionary<string, LinkedList<double>> ForComputingAverages = new Dictionary<string, LinkedList<double>>();
                foreach (Atom atom in AtomsInResidue)
                {

                    Atom[] CurrentBond = new Atom[2];
                    CurrentBond[0] = atom;
                    double[] ForNewOrigin = atom.GetRelativeAtomCoords();
                    foreach (Atom edge in atom.GetEdges())
                    {
                        double[] PosRelativeToRoot = { edge.GetRelativeAtomCoords()[0] - atom.GetRelativeAtomCoords()[0], edge.GetRelativeAtomCoords()[1] - atom.GetRelativeAtomCoords()[1], edge.GetRelativeAtomCoords()[2] - atom.GetRelativeAtomCoords()[2] };
                        Console.WriteLine("[\"" + atom.GetAtomName() + "" + edge.GetAtomName() + "\"] = new double[] { " + string.Join(",", PosRelativeToRoot) + " },");
                        CurrentBond[1] = edge;
                        ForComputingAverages = ComputeAverages(CurrentBond, ForComputingAverages);
                    }
                    // Console.Write("\n");
                }

            }
        }

        public double[] GetCartesian(double[] Spherical)
        {
            double x = Spherical[0] * Math.Sin(Spherical[1]) * Math.Cos(Spherical[2]);
            double y = Spherical[0] * Math.Sin(Spherical[1]) * Math.Sin(Spherical[2]);
            double z = Spherical[0] * Math.Cos(Spherical[1]);
            return new double[] { x, y, z };
        }


        public void PrintList()
        {

            //Console.WriteLine("I am trying to do this");

            foreach (string key in KeyValuePairsAtoms.Keys.ToList())
            {
                Console.Write(key + ": ");
                foreach (Atom edge in KeyValuePairsAtoms[key].GetEdges())
                {
                    Console.Write(edge.GetAtomName() + " || ");
                }

                Console.Write("\n");
            }
        }



        public double[] GetBondWeight(string BondPair)
        {
            return BondWeights[BondPair];
        }

    }
}
