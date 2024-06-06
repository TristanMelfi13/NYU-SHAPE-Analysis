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








    }
}
