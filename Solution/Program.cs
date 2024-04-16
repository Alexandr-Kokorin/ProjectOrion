using System;
using System.Linq.Dynamic.Core;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Linq.Expressions;
using System.IO;

namespace Solution
{
    internal class Program
    {

        static void Main(string[] args)
        {
            NelderMeadMethod method = new NelderMeadMethod("Math.Pow(x[0], 2) + x[0] * x[1] + Math.Pow(x[1], 2) - 6*x[0] -9*x[1]", 2);
            var result = method.Run();
            Console.WriteLine("Результат:");
            Console.WriteLine("Аргументы:");
            for (int i = 0; i < 2; i++) {
                Console.WriteLine(result[i]);
            }
            Console.WriteLine("Решение:");
            Console.WriteLine(method.Function(result));
            Console.Read();
        }
    }
}
