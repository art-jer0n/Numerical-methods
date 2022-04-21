using System;
namespace Two_dimensionalEquationOffThermalСonductivity
{
    class Program
    {
        static void Main(string[] args)
        {
            Solver solver = new Solver(100, 100, 100);
            solver.Run();
            Console.ReadKey();
        }

    }
}
