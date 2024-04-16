using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Dynamic.Core;
using System.Linq.Expressions;
using System.Text;
using System.Threading.Tasks;

namespace Solution
{
    internal class NelderMeadMethod
    {
        private const double edgeLengthInitial = 0.4; // Начальная длина ребра симплекса
        private const double edgeLengthLimit = 1.0e-5; // Предельное значение длины ребра симплекса
        private const double reflection = 1.0; // Коэффициент отражения симплекса
        private const double alpha = 2.0; // Коэффициент растяжения симплекса
        private const double beta = 0.5; // Коэффициент сжатия симплекса
        private const double gamma = 0.5; // Коэффициент редукции симплекса

        private readonly string expression; //"Math.Pow(x[0], 2) + x[0] * x[1] + Math.Pow(x[1], 2) - 6*x[0] -9*x[1]"
        private readonly double[] X; // Первая вершина начального симплекса (начальная точка)
        private readonly int count; // count - число аргументов функции

        private readonly double[,] simplex; // count + 1 - число вершин симплекса
        private readonly double[] functionValues;

        public NelderMeadMethod(string expression, int count) {
            this.expression = expression;
            this.count = count;
            X = new double[count];
            for (int i = 0; i < count; i++) {
                X[i] = 0;
            }
            simplex = new double[count, count + 1];
            functionValues = new double[count + 1];
        }

        // Выполняет поиск экстремума (минимума) функции F
        public double[] Run()
        {
            int imi = -1, ima = -1;
            int j = 0, kr = 0, jMx = 10000; // Предельное число шагов алгоритма (убрать после отладки)
            double[] X_R = new double[count];
            double Fmi, Fma, F_R;
            const int kr_todo = 60; // kr_todo - число шагов алгоритма, после выполнения которых симплекс восстанавливается

            MakeSimplex(X, edgeLengthInitial);
            while (IsStop() && j < jMx) {
                j++; // Число итераций
                kr++;
                if (kr == kr_todo) {
                    kr = 0;
                    SimplexRestore(); // Восстановление симплекса
                }
                Fmi = MinValue(ref imi);
                Fma = MaxValue(ref ima); // ima - Номер отражаемой вершины
                for (int i = 0; i < count; i++) {
                    X[i] = simplex[i, ima];
                }
                Reflection(ima); // Отражение
                for (int i = 0; i < count; i++) {
                    X_R[i] = simplex[i, ima];
                }
                F_R = Function(X_R); // Значение функции в вершине ima симплекса после отражения
                if (F_R > Fma) {
                    Compression(ima, Fma);
                }
                else if (F_R < Fmi) {
                    Stretching(ima, Fmi, F_R, X_R);
                }
                else {
                    functionValues[ima] = F_R;
                }
            }
            Console.WriteLine("Число итераций: " + j);
            return X;
        }

        private void Compression(int ima, double Fma)
        {
            double[] X2 = new double[count];
            double F_S;

            ShrinkingExpansion(ima, beta); // Сжатие
            for (int i = 0; i < count; i++) {
                X2[i] = simplex[i, ima];
            }
            F_S = Function(X2); // Значение функции в вершине ima симплекса после его сжатия
            if (F_S > Fma) {
                Reduction(ima, X2);
            }
            else {
                functionValues[ima] = F_S;
            }
        }

        private void Reduction(int ima, double[] X2)
        {
            for (int i = 0; i < count; i++) {
                simplex[i, ima] = X[i];
            }
            Reduction(ima); // Редукция
            for (int i = 0; i < count + 1; i++) {
                if (i == ima) continue;
                for (int j = 0; j < count; j++) {
                    X2[j] = simplex[j, i];
                }
                // Значения функций в вершинах симплекса после редукции. В вершине ima значение функции сохраняется
                functionValues[i] = Function(X2);
            }
        }

        private void Stretching(int ima, double Fmi, double F_R, double[] X_R)
        {
            double[] X2 = new double[count];
            double F_E;

            ShrinkingExpansion(ima, alpha); // Растяжение
            for (int i = 0; i < count; i++) {
                X2[i] = simplex[i, ima];
            }
            F_E = Function(X2); // Значение функции в вершине ima симплекса после его растяжения
            if (F_E > Fmi) {
                for (int j = 0; j < count; j++) {
                    simplex[j, ima] = X_R[j];
                }
                functionValues[ima] = F_R;
            }
            else {
                functionValues[ima] = F_E;
            }
        }

        public double Function(double[] x)
        {
            return (double)DynamicExpressionParser.ParseLambda(new[] { Expression.Parameter(typeof(double[]), "x") }, null, expression).Compile().DynamicInvoke(x);
        }

        // Создает из точки X регулярный симплекс с длиной ребра L и с NP + 1 вершиной
        // Формирует массив FN значений оптимизируемой функции F в вершинах симплекса
        private void MakeSimplex(double[] X, double L)
        {
            double qn = Math.Sqrt(1.0 + count) - 1.0;
            double q2 = L / Math.Sqrt(2.0) * count;
            double r1 = q2 * (qn + count);
            double r2 = q2 * qn;

            for (int i = 0; i < count; i++) {
                simplex[i, 0] = X[i];
            }
            for (int i = 1; i < count + 1; i++) {
                for (int j = 0; j < count; j++) {
                    simplex[j, i] = X[j] + r2;
                }
            }
            for (int i = 1; i < count + 1; i++) {
                simplex[i - 1, i] = simplex[i - 1, i] - r2 + r1;
            }
            for (int i = 0; i < count + 1; i++) {
                for (int j = 0; j < count; j++) {
                    X[j] = simplex[j, i];
                }
                functionValues[i] = Function(X); // Значения функции в вершинах начального симплекса
            }
        }

        private double[] CenterOfGravity(int k) // Центр тяжести симплекса
        {
            double s;
            double[] xc = new double[count];
            for (int i = 0; i < count; i++) {
                s = 0;
                for (int j = 0; j < count + 1; j++) {
                    s += simplex[i, j];
                }
                xc[i] = s;
            }
            for (int i = 0; i < count; i++) {
                xc[i] = (xc[i] - simplex[i, k]) / count;
            }
            return xc;
        }

        private void Reflection(int k) // Отражение вершины с номером k относительно центра тяжести
        {
            double[] xc = CenterOfGravity(k);
            for (int i = 0; i < count; i++) {
                simplex[i, k] = (1.0 + reflection) * xc[i] - simplex[i, k];
            }
        }

        private void Reduction(int k) // Редукция симплекса к вершине k
        {
            double[] xk = new double[count];
            for (int i = 0; i < count; i++) {
                xk[i] = simplex[i, k];
            }
            for (int j = 0; j < count; j++) {
                for (int i = 0; i < count; i++) {
                    simplex[i, j] = xk[i] + gamma * (simplex[i, j] - xk[i]);
                }
            }
            for (int i = 0; i < count; i++) {
                simplex[i, k] = xk[i]; // Восстанавливаем симплекс в вершине k
            }
        }

        private void ShrinkingExpansion(int k, double alpha_beta) // Сжатие/растяжение симплекса. alpha_beta – коэффициент растяжения/сжатия
        {
            double[] xc = CenterOfGravity(k);
            for (int i = 0; i < count; i++) {
                simplex[i, k] = xc[i] + alpha_beta * (simplex[i, k] - xc[i]);
            }
        }

        private double FindEdgeLength(double[] X2) // Длиина ребра симплекса
        {
            double L = 0;
            for (int i = 0; i < count; i++) {
                L += X2[i] * X2[i];
            }
            return Math.Sqrt(L);
        }

        private double MinValue(ref int imi) // Минимальный элемент массива и его индекс
        {
            double fmi = double.MaxValue;
            for (int i = 0; i < count + 1; i++) {
                if (functionValues[i] < fmi) {
                    fmi = functionValues[i];
                    imi = i;
                }
            }
            return fmi;
        }

        private double MaxValue(ref int ima) // Максимальный элемент массива и его индекс
        {
            double fma = double.MinValue;
            for (int i = 0; i < count + 1; i++) {
                if (functionValues[i] > fma) {
                    fma = functionValues[i];
                    ima = i;
                }
            }
            return fma;
        }

        private void SimplexRestore() // Восстанавление симплекса
        {
            int imi = -1, imi2 = -1;
            double fmi, fmi2 = double.MaxValue;
            double[] X = new double[count];
            double[] X2 = new double[count];
            fmi = MinValue(ref imi);

            for (int i = 0; i < count + 1; i++) {
                if (functionValues[i] != fmi && functionValues[i] < fmi2) {
                    fmi2 = functionValues[i];
                    imi2 = i;
                }
            }
            for (int i = 0; i < count; i++) {
                X[i] = simplex[i, imi];
                X2[i] = simplex[i, imi] - simplex[i, imi2];
            }
            MakeSimplex(X, FindEdgeLength(X2));
        }

        private bool IsStop() // Возвращает true, если длина хотя бы одного ребра симплекса превышает edgeLengthLimit, или false - в противном случае
        {
            double[] X = new double[count];
            double[] X2 = new double[count];
            for (int i = 0; i < count; i++) {
                for (int j = 0; j < count; j++) {
                    X[j] = simplex[j, i];
                }
                for (int j = i + 1; j < count + 1; j++) {
                    for (int k = 0; k < count; k++) {
                        X2[k] = X[k] - simplex[k, j];
                    }
                    if (FindEdgeLength(X2) > edgeLengthLimit) return true;
                }
            }
            return false;
        }
    }
}
