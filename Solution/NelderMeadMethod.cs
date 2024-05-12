using System;
using System.Linq.Dynamic.Core;
using System.Linq.Expressions;

namespace Solution
{
    public class NelderMeadMethod
    {
        private const double edgeLengthInitial = 0.4; // Начальная длина ребра симплекса
        private const double edgeLengthLimit = 1.0e-5; // Предельное значение длины ребра симплекса
        private const double reflection = 1.0; // Коэффициент отражения симплекса
        private const double alpha = 2.0; // Коэффициент растяжения симплекса
        private const double beta = 0.5; // Коэффициент сжатия симплекса
        private const double gamma = 0.5; // Коэффициент редукции симплекса

        public string expression; //"Math.Pow(x[1], 2) + x[1] * x[2] + Math.Pow(x[2], 2) - 6*x[1] -9*x[2]"
        public double[] X; // Первая вершина начального симплекса (начальная точка)
        public int count; // count - число аргументов функции

        public double[,] simplex; // count + 1 - число вершин симплекса
        public double[] functionValues;

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
            int imin = -1, imax = -1;
            int j = 0, kr = 0;
            double[] X_R = new double[count];
            double Fmin, Fmax, F_R;
            const int kr_todo = 60; // kr_todo - число шагов алгоритма, после выполнения которых симплекс восстанавливается

            MakeSimplex(X, edgeLengthInitial);
            while (IsStop()) {
                j++; // Число итераций
                kr++;
                if (kr == kr_todo) {
                    kr = 0;
                    SimplexRestore(); // Восстановление симплекса
                }
                Fmin = MinValue(ref imin);
                Fmax = MaxValue(ref imax); // imax - Номер отражаемой вершины
                for (int i = 0; i < count; i++) {
                    X[i] = simplex[i, imax];
                }
                Reflection(imax); // Отражение
                for (int i = 0; i < count; i++) {
                    X_R[i] = simplex[i, imax];
                }
                F_R = Function(X_R); // Значение функции в вершине imax симплекса после отражения
                if (F_R > Fmax) {
                    Compression(imax, Fmax);
                }
                else if (F_R < Fmin) {
                    Stretching(imax, Fmin, F_R, X_R);
                }
                else {
                    functionValues[imax] = F_R;
                }
            }
            Console.WriteLine("Число итераций: " + j);
            return X;
        }

        public void Compression(int imax, double Fmax)
        {
            double[] X2 = new double[count];
            double F_S;

            ShrinkingExpansion(imax, beta); // Сжатие
            for (int i = 0; i < count; i++) {
                X2[i] = simplex[i, imax];
            }
            F_S = Function(X2); // Значение функции в вершине imax симплекса после его сжатия
            if (F_S > Fmax) {
                Reduction(imax, X2);
            }
            else {
                functionValues[imax] = F_S;
            }
        }

        public void Reduction(int imax, double[] X2)
        {
            for (int i = 0; i < count; i++) {
                simplex[i, imax] = X[i];
            }
            Reduction(imax); // Редукция
            for (int i = 0; i < count + 1; i++) {
                if (i == imax) continue;
                for (int j = 0; j < count; j++) {
                    X2[j] = simplex[j, i];
                }
                // Значения функций в вершинах симплекса после редукции. В вершине imax значение функции сохраняется
                functionValues[i] = Function(X2);
            }
        }

        public void Stretching(int imax, double Fmin, double F_R, double[] X_R)
        {
            double[] X2 = new double[count];
            double F_E;

            ShrinkingExpansion(imax, alpha); // Растяжение
            for (int i = 0; i < count; i++) {
                X2[i] = simplex[i, imax];
            }
            F_E = Function(X2); // Значение функции в вершине imax симплекса после его растяжения
            if (F_E > Fmin) {
                for (int j = 0; j < count; j++) {
                    simplex[j, imax] = X_R[j];
                }
                functionValues[imax] = F_R;
            }
            else {
                functionValues[imax] = F_E;
            }
        }

        public double Function(double[] X)
        {
            double[] x = new double[count + 1];
            for (int i = 1; i < count + 1; i++) {
                x[i] = X[i - 1];
            }
            return (double)DynamicExpressionParser.ParseLambda(new[] { Expression.Parameter(typeof(double[]), "x") }, null, expression).Compile().DynamicInvoke(x);
        }

        // Создает из точки X регулярный симплекс с длиной ребра L и с NP + 1 вершиной
        // Формирует массив functionValues значений оптимизируемой функции F в вершинах симплекса
        public void MakeSimplex(double[] X, double L)
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

        public double[] CenterOfGravity(int k) // Центр тяжести симплекса
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

        public void Reflection(int k) // Отражение вершины с номером k относительно центра тяжести
        {
            double[] xc = CenterOfGravity(k);
            for (int i = 0; i < count; i++) {
                simplex[i, k] = (1.0 + reflection) * xc[i] - simplex[i, k];
            }
        }

        public void Reduction(int k) // Редукция симплекса к вершине k
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

        public void ShrinkingExpansion(int k, double alpha_beta) // Сжатие/растяжение симплекса. alpha_beta – коэффициент растяжения/сжатия
        {
            double[] xc = CenterOfGravity(k);
            for (int i = 0; i < count; i++) {
                simplex[i, k] = xc[i] + alpha_beta * (simplex[i, k] - xc[i]);
            }
        }

        public double FindEdgeLength(double[] X2) // Длиина ребра симплекса
        {
            double L = 0;
            for (int i = 0; i < count; i++) {
                L += X2[i] * X2[i];
            }
            return Math.Sqrt(L);
        }

        public double MinValue(ref int imin) // Минимальный элемент массива и его индекс
        {
            double fmin = double.MaxValue;
            for (int i = 0; i < count + 1; i++) {
                if (functionValues[i] < fmin) {
                    fmin = functionValues[i];
                    imin = i;
                }
            }
            return fmin;
        }

        public double MaxValue(ref int imax) // Максимальный элемент массива и его индекс
        {
            double fmax = double.MinValue;
            for (int i = 0; i < count + 1; i++) {
                if (functionValues[i] > fmax) {
                    fmax = functionValues[i];
                    imax = i;
                }
            }
            return fmax;
        }

        public void SimplexRestore() // Восстанавление симплекса
        {
            int imin = -1, imin2 = -1;
            double fmin, fmin2 = double.MaxValue;
            double[] X = new double[count];
            double[] X2 = new double[count];
            fmin = MinValue(ref imin);

            for (int i = 0; i < count + 1; i++) {
                if (functionValues[i] != fmin && functionValues[i] < fmin2) {
                    fmin2 = functionValues[i];
                    imin2 = i;
                }
            }
            for (int i = 0; i < count; i++) {
                X[i] = simplex[i, imin];
                X2[i] = simplex[i, imin] - simplex[i, imin2];
            }
            MakeSimplex(X, FindEdgeLength(X2));
        }

        public bool IsStop() // Возвращает true, если длина хотя бы одного ребра симплекса превышает edgeLengthLimit, или false - в противном случае
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
