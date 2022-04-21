using System;
namespace Quasi_LinearEquationOfThermalСonductivity
{
    class Program
    {
        static void Main(string[] args)
        {
            var solver = new Solver(10, 100);
            solver.Run();
            Console.ReadKey();
        }
    }
}

