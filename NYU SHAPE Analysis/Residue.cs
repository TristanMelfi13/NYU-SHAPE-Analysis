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
        
        public Residue(string residuename)
        {
            ResidueName = residuename;
        }

        public string GetName()
        {
            return ResidueName;
        }

        public LinkedList<Atom> GetAtomList()
        {
            return AtomsInResidue;
        }

        public void AddAtom(string[] AtomicInformation)
        {
            string CurrentRootAtom = AtomicInformation[1];
            string CurrentNeighbor = AtomicInformation[5];
            Atom atomRoot = null;
            Atom atomNeighbor = null;
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
                Atom NewRoot = new Atom(CurrentRootAtom, new double[] { double.Parse(AtomicInformation[2]), double.Parse(AtomicInformation[3]), double.Parse(AtomicInformation[4]) });
                Atom NewNeighbor = new Atom(CurrentNeighbor, new double[] { double.Parse(AtomicInformation[6]), double.Parse(AtomicInformation[7]), double.Parse(AtomicInformation[8]) });
                NewRoot.AddEdge(NewNeighbor);
                AtomsInResidue.AddLast(NewRoot);
                AtomsInResidue.AddLast(NewNeighbor);
                return;



            }
            else if (atomRoot != null && atomNeighbor == null)
            {
                Atom NewNeighbor = new Atom(CurrentNeighbor, new double[] { double.Parse(AtomicInformation[6]), double.Parse(AtomicInformation[7]), double.Parse(AtomicInformation[8]) });
                atomRoot.AddEdge(NewNeighbor);
                AtomsInResidue.AddLast(NewNeighbor);
                return;

            }
            else if (atomRoot == null && atomNeighbor != null)
            {
                Atom NewRoot = new Atom(CurrentRootAtom, new double[] { double.Parse(AtomicInformation[2]), double.Parse(AtomicInformation[3]), double.Parse(AtomicInformation[4]) });
                NewRoot.AddEdge(atomNeighbor);
                AtomsInResidue.AddLast(NewRoot);
                return;
            }



        }

        public void PrintList()
        {
            foreach(Atom atom in AtomsInResidue)
            {
                Console.Write(atom.GetAtomName() + " ---> ");
                foreach (Atom edge in atom.GetEdges())
                {
                    Console.Write(edge.GetAtomName() + " ");
                }
                Console.Write("\n");
            }
        }

    }
}
