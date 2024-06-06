using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class ExperimentalData
    {
        private string ReagentName = "";

        private Dictionary<string, LinkedList<Tuple<string, double>>> Results = new Dictionary<string, LinkedList<Tuple<string, double>>>();

        private Dictionary<string, double> AvgResults = new Dictionary<string, double>();

        public ExperimentalData(string name)
        {
            ReagentName = name;
        }

        public string GetReagentName()
        {
            return ReagentName;
        }

        public void AddToResults(string ReactionConcentration, LinkedList<Tuple<string, double>> ToAdd)
        {
            Results.Add(ReactionConcentration, ToAdd);
        }

        public Dictionary<string, LinkedList<Tuple<string, double>>> GetKeyValuePairs()
        {
            return Results;
        }

        public bool Check(string KeyToCheck)
        {
            return Results.ContainsKey(KeyToCheck);
        }


        public void CreateAvgDatas()
        {

            Func<double[], int, double> Avg = (v1, v2) => {
                double RunningAvg = 0;
                foreach(double d in v1)
                {
                    RunningAvg += d;
                }
                return RunningAvg / v2;
            };
            Dictionary<string, LinkedList<double>> ForCreatingAverages = new Dictionary<string, LinkedList<double>>();
            foreach (var kvp in Results)
            {
                LinkedList<Tuple<string, double>> Results = kvp.Value;
                // Console.WriteLine(string.Join("\n", Results));
                

                foreach (Tuple<string, double> kvp2 in Results)
                {
                    string NucleotideName = kvp2.Item1;


                    if (ForCreatingAverages.ContainsKey(NucleotideName))
                    {
                        ForCreatingAverages[NucleotideName].AddLast(kvp2.Item2);
                    } else
                    {
                        ForCreatingAverages.Add(NucleotideName, new LinkedList<double> { });
                        ForCreatingAverages[NucleotideName].AddLast(kvp2.Item2);
                    }
                }

            }
            foreach (var kvp in ForCreatingAverages)
            {
                double AverageReactivity = Avg(kvp.Value.ToArray(), kvp.Value.Count());
                Tuple<string, double> NucleotideValuePair = new Tuple<string, double>(kvp.Key, AverageReactivity);
                AvgResults.Add(kvp.Key, AverageReactivity);

            }


        }
    
        public Dictionary<string, double> GetAvgVals()
        {
            return AvgResults;
        }
    
    }
}
