using System;

namespace Solution
{
    internal class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Введите функцию: (Пример: x1^2+x1*x2+x2^2-6*x1-9*x2)");
            string input = Console.ReadLine();
            Parser parser = new Parser();
            int count = parser.getCount(input);
            NelderMeadMethod method = new NelderMeadMethod(parser.parse(input), count);
            var result = method.Run();
            Console.WriteLine("Результат:");
            Console.WriteLine(" Аргументы:");
            for (int i = 0; i < count; i++) {
                Console.WriteLine("  x" + (i+1) + " = " + result[i]);
            }
            Console.WriteLine(" Решение:");
            Console.WriteLine("  " + method.Function(result));
            Console.Read();
        }
    }
}
