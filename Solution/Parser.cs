using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace Solution
{
    internal class Parser
    {
        public string parse(string input)
        {
            string result = input;
            result = Regex.Replace(result, @"x(\d*)", "x[$1]");
            result = Regex.Replace(result, @"(x\[\d*\])\^(\d*)", "Math.Pow($1, $2)");
            result = Regex.Replace(result, @"\(([x+\-*/\d\[\]()]*])\)\^(\d*)", "Math.Pow($1, $2)");
            result = Regex.Replace(result, @"sin\((.*)\)", "Math.Sin($1)");
            result = Regex.Replace(result, @"cos\((.*)\)", "Math.Cos($1)");
            return result;
        }

        public int getCount(string input)
        {
            HashSet<string> result = new HashSet<string>();
            var matches = Regex.Matches(input, @"x\d*");
            foreach (var match in matches ) {
                result.Add(match.ToString());
            }
            return result.Count;
        }
    }
}
