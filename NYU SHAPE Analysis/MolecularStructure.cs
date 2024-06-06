using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class MolecularStructure
    {
        private string Name = "";

        private List<string> UnreactiveNucleotides = new List<string>();

        private List<string> ReactiveNucleotides = new List<string>();

        private List<string> HairpinNucleotides = new List<string>();

        private List<Tuple<string, string>> WCPairs = new List<Tuple<string, string>>();


        public MolecularStructure(string n)
        {
            Name = n;
        }

        public void AddToUnreactive(List<string> ToAdd)
        {
            lock (UnreactiveNucleotides)
            {
                var Current = UnreactiveNucleotides;
                UnreactiveNucleotides = Current.Union(ToAdd).ToList();
            }
        }

        public void AddToReactive(List<string> ToAdd)
        {
            lock (ReactiveNucleotides)
            {
                var Current = ReactiveNucleotides;
                ReactiveNucleotides = Current.Union(ToAdd).ToList();
            }
        }

        public void AddToHairpin(List<string> ToAdd)
        {
            lock (HairpinNucleotides)
            {
                var Current = HairpinNucleotides;
                HairpinNucleotides = Current.Union(ToAdd).ToList();
            }
        }
    }
}
