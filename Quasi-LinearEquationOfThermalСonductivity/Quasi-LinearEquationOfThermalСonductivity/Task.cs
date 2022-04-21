using System;
namespace Quasi_LinearEquationOfThermalСonductivity
{
    class Task
    {
        public static double U(double x, double t)
        {
            return Math.Pow(x, 2) * Math.Pow(t, 2) + x + 1;
        }
        public static double Fi(double u)
        {
            return Math.Pow(u, 3);
        }
        public static double DFi(double u)
        {
            return 3 * Math.Pow(u, 2);
        }
        public static double BoundaryCondition_m1()
        {
            return 1.0;
        }
        public static double BoundaryCondition_m2(double t)
        {
            return Math.Pow(t, 2) + 2.0;
        }
        public static double BoundaryCondition_u0(double x)
        {
            return x + 1.0;
        }
        public static double F(double x, double t)
        {
            return 6 * t * Math.Pow(x, 2) * Math.Pow((Math.Pow(t, 2) * Math.Pow(x, 2) + x + 1), 2) - 2 * Math.Pow(t, 2);
        }
        public static double MaximumElem(double[,] array, double N, double M)
        {
            double result = 0.0;
            for (int i = 0; i < N + 1; i++)
            {
                for (int j = 0; j < M + 1; j++)
                    if (array[i, j] > result) result = array[i, j];
            }
            return result;

        }
        public static void Error(double[,] y, double h, double tay, double[,] z, double N, double M)
        {
            for (int i = 0; i < N + 1; i++)
                for (int j = 0; j < M + 1; j++)
                    z[i, j] = Math.Abs((y[i, j] - U(h * i, tay * j)) / y[i, j]);
        }
        public static double Norm(double[,] y1, double[,] y2, int j, double N)
        {
            double result = 0.0;
            for (int i = 0; i < N + 1; i++)
            {
                double z = Math.Abs(y2[i, j] - y1[i, j]);
                if (z > result) result = z;
            }
            return result;
        }
    }

}
