using System;
namespace Quasi_LinearEquationOfThermalСonductivity
{
    class Solver
    {
        public double eps = 10e-4;
        public int N, M;
        public double h, tay;
        public double[,] y, z, a;

        public Solver(int N, int M)
        {
            this.N = N;
            this.M = M;
            h = 1.0 / N;
            tay = 1.0 / M;
        }

        private void InitialConditions()
        {
            y = new double[N + 1, M + 1];
            a = new double[N + 1, M + 1];
            z = new double[N + 1, M + 1];

            for (int i = 0; i <= N; i++)
                y[i, 0] = Task.BoundaryCondition_u0(h * i);

            for (int j = 0; j <= M; j++)
            {
                y[0, j] = Task.BoundaryCondition_m1();
                y[N, j] = Task.BoundaryCondition_m2(tay * j);
            }
            for (int i = 0; i < N + 1; i++)
                for (int j = 0; j < M + 1; j++)
                    a[i, j] = y[i, j];
        }
        public void Run()
        {
            InitialConditions();
            double[] alpha = new double[N + 1];
            double[] beta = new double[N + 1];

            for (int j = 0; j < M; j++)
            {
                alpha[1] = 0.0;
                beta[1] = Task.BoundaryCondition_m1();
                do
                {
                    for (int k = 0; k < N + 1; k++) y[k, j + 1] = a[k, j + 1];
                    for (int i = 1; i < N; i++)
                    {
                        double A = tay / h / h;
                        double C = 2.0 * tay / h / h + Task.DFi(y[i, j + 1]);
                        double B = tay / h / h;
                        double F = Task.Fi(a[i, j]) - Task.Fi(y[i, j + 1]) + Task.DFi(y[i, j + 1]) * y[i, j + 1] + tay * Task.F(h * i, tay * j);

                        alpha[i + 1] = B / (C - A * alpha[i]);
                        beta[i + 1] = (A * beta[i] + F) / (C - A * alpha[i]);
                    }
                    for (int i = N - 1; i > 0; i--)
                        a[i, j + 1] = alpha[i + 1] * a[i + 1, j + 1] + beta[i + 1];
                }
                while (Task.Norm(y, a, j + 1, N) >= eps);
            }
            Task.Error(a, h, tay, z, N, M);
            Console.WriteLine("\n N:{0}\t M:{1}\t max. error: {2}", N, M, (Task.MaximumElem(z, N, M)));
        }
    }
}
