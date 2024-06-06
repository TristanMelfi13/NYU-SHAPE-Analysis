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

        private string NucleotideWeBelongTo = "";


        private string VisitedColor = "";

        private double[] XYZCoordinates = new double[3];

        private double[] RelativeCoordinates = new double[3]; // Really just the distance between the other bonds

        private Dictionary<Atom, double[]> EdgeWeightsDistance = new Dictionary<Atom, double[]>();

        private LinkedList<Atom> AtomicEdges = new LinkedList<Atom>();

        public Atom(string name, double[] Coords, string nucleotideWeBelongTo)
        {
            AtomName = name;
            RelativeCoordinates = Coords; //new double[] { 0.0, 0.0, 0.0 };
            NucleotideWeBelongTo = nucleotideWeBelongTo;
        }

        public string GetNucleotideName()
        {
            return NucleotideWeBelongTo;
        }


        public double[] GetAtomCoords()
        {
            return RelativeCoordinates;
        }

        public void AddEdge(Atom ToAdd)
        {
            AtomicEdges.AddLast(ToAdd);
        }



        public void AddEdge(Atom ToAdd, double[] Weight)
        {
            AtomicEdges.AddLast(ToAdd);

            EdgeWeightsDistance.Add(ToAdd, Weight);



        }

        public LinkedList<Atom> GetEdges()
        {
            List<Atom> ToReturn = AtomicEdges.Distinct().ToList();
            return new LinkedList<Atom>(ToReturn);
        }


        public string GetAtomName() { return AtomName; }


        public double[] GetRelativeAtomCoords() { return RelativeCoordinates; }
    
    }
}
