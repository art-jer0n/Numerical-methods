using System;
using System.Collections.Generic;
using System.Text;

namespace Two_dimensionalEquationOffThermalСonductivity
{
    class Solver
    {
        private int N1, N2, M;
        private double h1, h2, tay;
        private double[,,] y, z;

        public Solver(int N1, int N2, int M)
        {
            this.N1 = N1;
            this.N2 = N2;
            this.M = M;

            h1 = 1.0 / N1;
            h2 = 1.0 / N2;
            tay = 1.0 / M;
        }

        private void InitialConditions()
        {
            y = new double[N1 + 1, N2 + 1, M + 1];
            z = new double[N1 + 1, N2 + 1, M + 1];

            for (int i = 0; i <= N1; i++)
                for (int j = 0; j <= N2; j++)
                    y[i, j, 0] = Task.u0(h1 * i, h2 * j);

            for (int j = 0; j <= N2; j++)
                for (int n = 1; n <= M; n++)
                {
                    y[0, j, n] = Task.m1(h2 * j, tay * n);
                    y[N1, j, n] = Task.m2(h2 * j, tay * n);
                }

            for (int i = 0; i <= N1; i++)
                for (int n = 1; n <= M; n++)
                {
                    y[i, 0, n] = Task.m3(h1 * i, tay * n);
                    y[i, N2, n] = Task.m4(h1 * i, tay * n);
                }
        }
        public void Run()
        {
            InitialConditions();
            for (int n = 0; n < M; n++)
            {
                double[,] a = new double[N1 + 1, N2 + 1];

                for (int j = 0; j <= N2; j++)
                {
                    a[0, j] = (Task.m1(h2 * j, tay * (n + 1)) + Task.m1(h2 * j, tay * n)) / 2.0
                        - tay / 4.0 * ((Task.m1(h2 * (j + 1), tay * (n + 1)) - 2.0 * Task.m1(h2 * j, tay * (n + 1)) + Task.m1(h2 * (j - 1), tay * (n + 1))) / Math.Pow(h1, 2)
                                       - (Task.m1(h2 * (j + 1), tay * n) - 2.0 * Task.m1(h2 * j, tay * n) + Task.m1(h2 * (j - 1), tay * n)) / Math.Pow(h1, 2));


                    a[N1, j] = (Task.m2(h2 * j, tay * (n + 1)) + Task.m2(h2 * j, tay * n)) / 2.0
                        - tay / 4.0 * ((Task.m2(h2 * (j + 1), tay * (n + 1)) - 2.0 * Task.m2(h2 * j, tay * (n + 1)) + Task.m2(h2 * (j - 1), tay * (n + 1))) / Math.Pow(h1, 2)
                                       - (Task.m2(h2 * (j + 1), tay * n) - 2.0 * Task.m2(h2 * j, tay * n) + Task.m2(h2 * (j - 1), tay * n)) / Math.Pow(h1, 2));
                }

                for (int j = 1; j < N2; j++)
                {
                    double[] alpha = new double[N1 + 1];
                    double[] beta = new double[N1 + 1];

                    alpha[1] = 0.0;
                    beta[1] = a[0, j];

                    for (int i = 1; i < N1; i++)
                    {
                        double A = 1.0 / Math.Pow(h1, 2);
                        double B = 1.0 / Math.Pow(h1, 2);
                        double C = 2.0 * (1.0 / Math.Pow(h1, 2) + 1.0 / tay);
                        double F = 2.0 / tay * y[i, j, n] + (y[i, j + 1, n] - 2.0 * y[i, j, n] + y[i, j - 1, n]) / Math.Pow(h2, 2) + Task.F(h1 * i, h2 * j, tay * (n + 1));
                        alpha[i + 1] = B / (C - A * alpha[i]);
                        beta[i + 1] = (A * beta[i] + F) / (C - A * alpha[i]);

                    }
                    for (int i = N1 - 1; i > 0; i--)
                        a[i, j] = alpha[i + 1] * a[i + 1, j] + beta[i + 1];

                }

                for (int i = 1; i < N1; i++)
                {
                    double[] alpha = new double[N2 + 1];
                    double[] beta = new double[N2 + 1];

                    alpha[1] = 0.0;
                    beta[1] = y[i, 0, n + 1];

                    for (int j = 1; j < N2; j++)
                    {
                        double A = 1.0 / Math.Pow(h2, 2);
                        double B = 1.0 / Math.Pow(h2, 2);
                        double C = 2.0 * (1.0 / Math.Pow(h2, 2) + 1.0 / tay);
                        double F = 2.0 / tay * a[i, j] + 1.0 / Math.Pow(h1, 2) * (a[i + 1, j] - 2.0 * a[i, j] + a[i - 1, j]) + Task.F(h1 * i, h2 * j, tay * n);

                        alpha[j + 1] = B / (C - A * alpha[j]);
                        beta[j + 1] = (A * beta[j] + F) / (C - A * alpha[j]);

                    }
                    for (int j = N2 - 1; j > 0; j--)
                        y[i, j, n + 1] = alpha[j + 1] * y[i, j + 1, n + 1] + beta[j + 1];

                }
            }
            Task.Error(y, h1, h2, tay, z, N1, N2, M);
            Console.WriteLine(" N1:{0,-5} N2:{1,-5} M:{2,-5} error: {3}", N1, N2, M, Task.Max(z, N1, N2, M));

        }
    }

}
