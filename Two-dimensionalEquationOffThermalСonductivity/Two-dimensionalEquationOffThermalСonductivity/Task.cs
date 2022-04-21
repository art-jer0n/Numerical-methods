using System;
using System.Collections.Generic;
using System.Text;

namespace Two_dimensionalEquationOffThermalСonductivity
{
    class Task
    {
        public static double U(double x1, double x2, double t)
        {
            return (Math.Pow(x1, 4) + Math.Pow(x2, 4)) * Math.Pow(t, 2) + x1 + x2;
        }
        public static double m1(double x2, double t)
        {
            return Math.Pow(x2, 4) * Math.Pow(t, 2) + x2;
        }
        public static double m2(double x2, double t)
        {
            return (Math.Pow(x2, 4) + 1.0) * Math.Pow(t, 2) + x2 + 1.0;
        }
        public static double m3(double x1, double t)
        {
            return Math.Pow(x1, 4) * Math.Pow(t, 2) + x1;
        }
        public static double m4(double x1, double t)
        {
            return (Math.Pow(x1, 4) + 1.0) * Math.Pow(t, 2) + x1 + 1.0;
        }
        public static double u0(double x1, double x2)
        {
            return x1 + x2;
        }
        public static double F(double x1, double x2, double t)
        {
            return 2.0 * t * (Math.Pow(x1, 4) + Math.Pow(x2, 4)) - 12.0 * Math.Pow(t, 2) * (Math.Pow(x1, 2) + Math.Pow(x2, 2));
        }
        public static double Max(double[,,] array, double N1, double N2, double M)
        {
            double result = 0.0;
            for (int i = 0; i < N1; i++)
                for (int j = 0; j < N2; j++)
                    for (int n = 0; n < M; n++)
                        if (array[i, j, n] > result) result = array[i, j, n];
            return result;
        }
        public static void Error(double[,,] array, double h1, double h2, double tay, double[,,] z, double N1, double N2, double M)
        {
            for (int i = 0; i < N1; i++)
                for (int j = 0; j < N2; j++)
                    for (int n = 0; n < M; n++)
                        z[i, j, n] = Math.Abs((array[i, j, n] - U(h1 * i, h2 * j, tay * n)) / array[i, j, n]);
        }
    }

}
