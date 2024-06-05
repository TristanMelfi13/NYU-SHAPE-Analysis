using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class Atom
    {
        private string AtomName = "";
        private double[] Coordinates = new double[3];
        private LinkedList<Atom> AtomicEdges = new LinkedList<Atom>();

        public Atom(string name, double[] Coords)
        {
            AtomName = name;
            Coordinates = Coords;
        }

        public void AddEdge(Atom ToAdd)
        {
            AtomicEdges.AddLast(ToAdd);
        }

        public LinkedList<Atom> GetEdges()
        {
            return AtomicEdges;
        }


        public string GetAtomName() { return AtomName; }

        public double[] GetAtomCoords() { return Coordinates; }
    }
}
