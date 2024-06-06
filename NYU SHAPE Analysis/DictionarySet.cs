using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NYU_SHAPE_Analysis
{
    internal class DictionarySet<T, K>
        where T : notnull
        where K : notnull
    {
        private string ListName = "";

        private Dictionary<string, Dictionary<T, List<K>>> Dictionary_Set = new Dictionary<string, Dictionary<T, List<K>>>();

        public DictionarySet(string listName) // Initilize
        {
            ListName = listName;
        }

        private void AddNewDictionary(string name)
        {
            lock (this)
            {
                Dictionary_Set.Add(name, new Dictionary<T, List<K>>());
            }
        }


        public void AddToDictionary(string DictionaryName, T Key, K Value)
        {
            lock (this)
            {

                if (!Dictionary_Set.ContainsKey(DictionaryName))
                {
                    AddNewDictionary(DictionaryName);
                }
                if (!Dictionary_Set[DictionaryName].ContainsKey(Key))
                {
                    Dictionary_Set[DictionaryName].Add(Key, new List<K>());
                }
                Dictionary_Set[DictionaryName][Key].Add(Value);
            }
        }


        public Dictionary<T, List<K>> GetSubSet(string name)
        {
            return Dictionary_Set[name];
        }

        public List<string> GetSuperKeys()
        {
            return Dictionary_Set.Keys.ToList();
        }

        public Dictionary<T, List<K>> GetDictionary(string name)
        {
            if (Dictionary_Set.ContainsKey(name))
            {
                return Dictionary_Set[name];
            }
            return null; 
        }

        public string GetSetName()
        {
            return ListName;
        }

        public Dictionary<string, Dictionary<T, List<K>>> GetSet()
        {
            return Dictionary_Set;
        }

        public bool ContainsSuperKey(string name) 
        {
            return Dictionary_Set.ContainsKey(name);
        }

        public void SetUpNewList(string SuperName)
        {
            Dictionary_Set.Add(SuperName, new Dictionary<T, List<K>>());
        }

        public bool ContainsKey(string SupName, T Name)
        {
            if (Dictionary_Set.ContainsKey(SupName) && Dictionary_Set[SupName].ContainsKey(Name))
            {
                return true;
            } else
            {
                return false;
            }
        }


        public List<string> ListDictionaryNames()
        {
            return new List<string>(Dictionary_Set.Keys);
        }

        public string GetName()
        {
            return ListName;
        }

        public void PrintSet()
        {
            foreach (var d in Dictionary_Set)
            {
                Console.WriteLine("=====" + d.Key + "=====");
                foreach (var kvp in d.Value)
                {
                    Console.WriteLine(kvp.Key + " ----> {" + string.Join(", ", kvp.Value) + "}");
                }
            }
        }

        public void Average()
        {
            foreach (var d in Dictionary_Set)
            {
                foreach (var kvp in d.Value)
                {
                    double RunningSum = 0;
                }
            }
        }

        

    }
}
