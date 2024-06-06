using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class ConsolidatedResults
    {
        private string Name = "";
        private DictionarySet<string, double> Results = new DictionarySet<string, double>("");


        public ConsolidatedResults(string name) 
        {
            Name = name;
        }

        public string GetName()
        {
            return Name;
        }

        public DictionarySet<string, double> GetDictionarySet()
        {
            return Results;
        }

        private string GetSelfNucleotideType(string input)
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

        private double Avg(List<double> ToAvg)
        {
            double RunningSum = 0;
            foreach (var d in ToAvg)
            {
                RunningSum += d;
            }
            return RunningSum / ToAvg.Count();
        }

        public void AddToResults(DictionarySet<string, double> ToConsolidate)
        {
            lock (Results) 
            {
                foreach (var Interaction in ToConsolidate.GetSet())
                {
                    var Split = Interaction.Key.Split('\t');
                    var Resi1 = GetSelfNucleotideType(Split[0]);
                    var Resi2 = GetSelfNucleotideType(Split[1]);
                    foreach(var kvp in Interaction.Value)
                    {
                        // Console.WriteLine(Resi1 + Resi2 + " " + kvp.Key + " " + string.Join(" ", kvp.Value));
                        Results.AddToDictionary(Resi1 + Resi2, kvp.Key, Avg(kvp.Value));
                    }
                }    
            }
        }


        public void PrintResults()
        {
            Results.PrintSet();
        }

    }
}
