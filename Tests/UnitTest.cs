using Microsoft.VisualStudio.TestTools.UnitTesting;
using Solution;
using System.Linq;

namespace Tests
{
    [TestClass]
    public class UnitTest
    {
        [TestMethod]
        public void TestMethod1()
        {
            /// Тест создания экземпляра класса
            var expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            var count = 2;
            var nm = new NelderMeadMethod(expression, count);
            //проверяем не является ли nm null - успешно ли создан экземпляр класса
            Assert.IsNotNull(nm);
        }

        [TestMethod]
        public void TestMethod2()
        {
            //Тест проверки функции Function
            NelderMeadMethod method = new NelderMeadMethod("Math.Pow(x[1], 2) + Math.Pow(x[2], 2)", 2);
            double[] inputs = new double[] { 1.0, 2.0 };
            double expected = 5.0;
            double actual = method.Function(inputs);
            //проверяем значение функции в точке
            Assert.AreEqual(expected, actual);
        }

        [TestMethod]
        public void TestMethod3()
        {
            //Тест проверки функции MakeSimplex
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.5;//длина ребра
            nmm.MakeSimplex(X, L);
            double[,] expectedSimplex = new double[count, count + 1];
            expectedSimplex[0, 0] = 1.0;
            expectedSimplex[1, 0] = 2.0;
            expectedSimplex[0, 1] = 2.9;
            expectedSimplex[1, 1] = 2.5;
            expectedSimplex[0, 2] = 1.5;
            expectedSimplex[1, 2] = 3.9;
            for (int i = 0; i < count+1; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    Assert.AreEqual(expectedSimplex[j, i], nmm.simplex[j, i], 0.1);
                }
            }
        }

        [TestMethod]
        public void TestMethod4()
        {
            //Тест проверки функции MaxValue
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2) + x[1] * x[2] - 6*x[1] -9*x[2]";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.4;
            //находим значение функции в каждой точке симплекса
            nmm.MakeSimplex(X, L);
            int ima = -1;
            double expected = -17;
            double actual = nmm.MaxValue(ref ima);
            Assert.AreEqual(expected, actual);
            Assert.AreEqual(0, ima);
        }

        [TestMethod]
        public void TestMethod5()
        {
            //Тест проверки функции MinValue
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2) + x[1] * x[2]";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.4;
            nmm.MakeSimplex(X, L);
            int imi = -1;
            double expected = 7;
            //находим значение функции в каждой точке симплекса
            double actual = nmm.MinValue(ref imi);
            Assert.AreEqual(expected, actual);
            Assert.AreEqual(0, imi);
        }

        [TestMethod]
        public void TestMethod6()
        {
            //Тест проверки функции 
            string expression = "Math.Pow(x[1], 2) + x[1] * x[2] + Math.Pow(x[2], 2) - 6*x[1] -9*x[2]";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            //double[] X = new double[] { 1.0, 2.0 };
            //double L = 0.4;
            //nmm.MakeSimplex(X, L);
            double[] X2 = new double[] { 1.2, 2.2 };
            double expected = 2.5;
            double actual = nmm.FindEdgeLength(X2);
            Assert.AreEqual(expected, actual, 0.01);
        }

        [TestMethod]
        public void TestMethod7()
        {
            //Тест на проверку функции СenterOfGravity
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.5;
            //создали симплекс
            nmm.MakeSimplex(X, L);
            int k = 0;
            double[] expected = new double[] { 2.2, 3.2 };
            double[] actual = nmm.CenterOfGravity(k);
            Assert.AreEqual(expected[0], actual[0], 0.1);
            Assert.AreEqual(expected[1], actual[1], 0.1);
        }

        [TestMethod]
        public void TestMethod8()
        {
            //Тест на проверку функции IsStop
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.5;
            nmm.MakeSimplex(X, L);
            bool expected = true;
            bool actual = nmm.IsStop();
            Assert.AreEqual(expected, actual);
        }

        [TestMethod]
        public void TestMethod9()
        {
            //Тест на проверку функции SimplexRestore
            string expression = "Math.Pow(x[1], 2) + x[1] * x[2] + Math.Pow(x[2], 2) - 6*x[1] -9*x[2]";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.4;
            nmm.MakeSimplex(X, L);
            nmm.SimplexRestore();
            double[] expected = new double[] { 1.4, 7.5 };
            double[] actual = nmm.simplex.Cast<double>().ToArray();
            Assert.AreEqual(expected[0], actual[0], 0.1);
            Assert.AreEqual(expected[1], actual[1], 0.1);
        }

        [TestMethod]
        public void TestMethod10()
        {
            //Тест на проверку функции ShrinkingExpansion
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.5;
            nmm.MakeSimplex(X, L);
            int k = 0;
            double alpha_beta = 2.0;
            double[] expected = new double[] { -0.2, 2.9 };
            nmm.ShrinkingExpansion(k, alpha_beta);
            double[] actual = nmm.simplex.Cast<double>().ToArray();
            Assert.AreEqual(expected[0], actual[0], 0.1);
            Assert.AreEqual(expected[1], actual[1], 0.1);
        }

        [TestMethod]
        public void TestMethod11()
        {
            //Тест на проверку функции Reduction
            string expression = "Math.Pow(x[1], 2) + x[1] * x[2] + Math.Pow(x[2], 2) - 6*x[1] -9*x[2]";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.4;
            nmm.MakeSimplex(X, L);
            int k = 0;
            nmm.Reduction(k);
            double[] expected = new double[] { 1.0, 1.7 };
            double[] actual = nmm.simplex.Cast<double>().ToArray();
            Assert.AreEqual(expected[0], actual[0], 0.1);
            Assert.AreEqual(expected[1], actual[1], 0.1);
        }

        [TestMethod]
        public void TestMethod12()
        {
            //Тест на проверку функции Reflection
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.4;
            nmm.MakeSimplex(X, L);
            int k = 0;
            double[] expected = new double[] { 2.9, 2.5 };
            nmm.Reflection(k);
            double[] actual = nmm.simplex.Cast<double>().ToArray();
            Assert.AreEqual(expected[0], actual[0], 0.1);
            Assert.AreEqual(expected[1], actual[1], 0.1);
        }

        [TestMethod]
        public void TestMethod13()
        {
            //Тест на проверку функции Compression
            string expression = "Math.Pow(x[1], 2) + Math.Pow(x[2], 2)";
            int count = 2;
            NelderMeadMethod nmm = new NelderMeadMethod(expression, count);
            double[] X = new double[] { 1.0, 2.0 };
            double L = 0.5;
            nmm.MakeSimplex(X, L);
            int ima = -1;
            double Fma = nmm.MaxValue(ref ima);
            nmm.Compression(ima, Fma);
            double[,] expectedSimplex = new double[count, count + 1];
            expectedSimplex[0, 0] = 1.0;
            expectedSimplex[1, 0] = 2.0;
            expectedSimplex[0, 1] = 2.9;
            expectedSimplex[1, 1] = 2.5;
            expectedSimplex[0, 2] = 1.7;
            expectedSimplex[1, 2] = 3.0;
            for (int i = 0; i < count+1; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    Assert.AreEqual(expectedSimplex[j, i], nmm.simplex[j, i], 0.1);
                }
            }
        }

        [TestMethod]
        public void TestMethod14()
        {
            // Тест на проверку функции Reduction
            int count = 2;
            NelderMeadMethod nm = new NelderMeadMethod("x[1] * x[1] + x[2] * x[2]", count);
            double[] X = new double[] { 1, 1 };
            nm.MakeSimplex(X, 0.4);
            int ima = 1;
            double[] X2 = new double[] { 2, 2 };
            nm.Reduction(ima, X2);
            double[,] expectedSimplex = new double[nm.count,nm.count+1];
            expectedSimplex[0, 0] = 0.5;
            expectedSimplex[1, 0] = 0.5;
            expectedSimplex[0, 1] = 0;
            expectedSimplex[1, 1] = 0;
            expectedSimplex[0, 2] = 1.4;
            expectedSimplex[1, 2] = 2.5;
            for (int i = 0; i < nm.count+1; i++)
            {
                for (int j = 0; j < nm.count; j++)
                {
                    Assert.AreEqual(expectedSimplex[j, i], nm.simplex[j, i], 0.1);
                }
            }
            double[] expectedFunctionValues = new double[nm.count + 1];
            expectedFunctionValues[0] = 0.5;
            expectedFunctionValues[1] = 8.4;
            expectedFunctionValues[2] = 8.4;
            for (int i = 0; i < nm.count + 1; i++)
            {
                Assert.AreEqual(expectedFunctionValues[i], nm.functionValues[i], 0.1);
            }
        }

        [TestMethod]
        public void TestMethod15()
        {
            //Тест на проверку функции Stretching
            NelderMeadMethod nm = new NelderMeadMethod("x[1] * x[1] + x[2] * x[2]", 2);
            double[] X = new double[] { 1, 1 };
            nm.MakeSimplex(X, 0.4);
            int ima = 1;
            double Fmi = 2;
            double F_R = 1;
            double[] X_R = new double[] { 2, 2 };
            nm.Stretching(ima, Fmi, F_R, X_R);
            for (int i = 0; i < 2; i++)
            {
                Assert.AreEqual(X_R[i], nm.simplex[i, ima]);
            }
            Assert.AreEqual(F_R, nm.functionValues[ima]);
        }

        [TestMethod]
        public void TestMethod16()
        {
            //Математический тест номер 1
            var nm = new NelderMeadMethod("(x[1] - 2) * (x[1] - 2) + (x[2] - 3) * (x[2] - 3)", 2);
            double[] X = new double[] { 1, 1 };
            double[] result = nm.Run();
            Assert.AreEqual(1.9, result[0], 0.1);
            Assert.AreEqual(2.9, result[1], 0.1);
            Assert.AreEqual(0, nm.Function(result), 0.1);
        }

        [TestMethod]
        public void TestMethod17()
        {
            //Математический тест номер 2
            var nm = new NelderMeadMethod("Math.Pow(x[1]-2, 2) + Math.Pow(x[2]-3, 2)", 2);
            double[] X = new double[] { 0, 0 };
            double[] result = nm.Run();
            Assert.AreEqual(1.9, result[0], 0.1);
            Assert.AreEqual(2.9, result[1], 0.1);
            Assert.AreEqual(0, nm.Function(result), 0.1);
        }

        [TestMethod]
        public void TestMethod18()
        {
            //Математический тест номер 3
            var nm = new NelderMeadMethod("Math.Pow(x[1], 3) + Math.Pow(x[2], 3) - 3*x[1]*x[2]", 2);
            double[] X = new double[] { 1, 1 };
            double[] result = nm.Run();
            Assert.AreEqual(1, result[0], 0.1);
            Assert.AreEqual(1, result[1], 0.1);
            Assert.AreEqual(-1, nm.Function(result), 0.1);
        }

        [TestMethod]
        public void TestMethod19()
        {
            //Математический тест номер 4
            var nm = new NelderMeadMethod("Math.Sin(x[1]) + Math.Cos(x[2])", 2);
            double[] X = new double[] { 0, 0 };
            double[] result = nm.Run();
            Assert.AreEqual(-1.5, result[0], 0.1);
            Assert.AreEqual(3.14, result[1], 0.1);
            Assert.AreEqual(-1.9, nm.Function(result), 0.1);
        }

        [TestMethod]
        public void TestMethod20()
        {
            //Математический тест номер 5
            var nm = new NelderMeadMethod("Math.Pow(x[1], 2) + Math.Pow(x[2], 2) - 10 * Math.Cos(x[1] * x[2])", 2);
            double[] X = new double[] { 1, 1 };
            double[] result = nm.Run();
            Assert.AreEqual(0, result[0], 0.1);
            Assert.AreEqual(0, result[1], 0.1);
            Assert.AreEqual(-10, nm.Function(result), 0.1);
        }

        [TestMethod]
        public void TestMethod21()
        {
            //Математический тест номер 6
            var method = new NelderMeadMethod("Math.Exp(Math.Pow(x[1], 2) + Math.Pow(x[2], 2))", 2);
            method.X = new double[] { 1, 1 };
            var result = method.Run();
            Assert.AreEqual(0, result[0], 0.1);
            Assert.AreEqual(0, result[1], 0.1);
            Assert.AreEqual(1, method.Function(result), 0.1);
        }
    }
}
