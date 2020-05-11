using System;

namespace Integral
{
    class Program
    {
        static double f(double x)
            => 2 * Math.Cos(2.5 * x) * Math.Exp(x * 1.0 / 3) + 4 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + x;
        static double NewtonKots(double lim1, double lim2, double step)
        {
            double integral = 0;
            while(lim2 <= b)
            {
                double[] nodes;
                double[] moments;
                moments = new double[3];

                // find moments
                moments[0] = functions.IntegralX0(lim1,lim2);
                moments[1] = functions.IntegralX1(lim1,lim2);
                moments[2] = functions.IntegralX2(lim1,lim2);

                // add 3 new nodes
                nodes = new double[3];
                nodes[0] = lim1;
                nodes[1] = lim1 + (lim2-lim1)/2;
                nodes[2] = b;
                
                // declare the coeff of system 
                double[,] matrixArray = new double[3,3];
                matrixArray[0,0] = 1;
                matrixArray[0,1] = 1;
                matrixArray[0,2] = 1;
                matrixArray[1,0] = nodes[0];
                matrixArray[1,1] = nodes[1];
                matrixArray[1,2] = nodes[2];
                matrixArray[2,0] = nodes[0] * nodes[0];
                matrixArray[2,1] = nodes[1] * nodes[1];
                matrixArray[2,2] = nodes[2] * nodes[2];

                // declare matrix
                SquareMatrix matrix = new SquareMatrix(matrixArray,3);
                matrix.LUDecomposition();

                // find the solution
                var solv = matrix.SolutionSystem(moments);

                for(int i = 0;i < 3;i++)
                    integral += solv[0,i] * f(nodes[i]);

                lim1 = lim2;
                lim2 += step;
            }
            return integral;
        }
        
        static double a = 1.5;
        static    double b = 3.3;
        static    double alpha = 1/3;
        static    double beta = 0;
        static    double exactValue = 7.077031437995793610263911711602477164432;
        static    Functions functions = new Functions();
        static void Main(string[] args)
        {

            // * * * * * * * * * * * * * Вариант ньютона - Котса * * * * * * * * * * * *
            Console.WriteLine("Вариант Ньютона-Котса : " + "\n");

            
            double integral = NewtonKots(a,b,1);
            double error = 0.6406;

            WriteText("Вариант Ньютона-Котса :",integral,Math.Abs(exactValue - integral),error);
            
            //* * * * * * * * * * * * * Cоставная ИКФ с нужной точностью* * * * * * * * * *
            double h = Math.Ceiling((b - a) / 10);
            double lim1 = a;
            double step = (b - a) / 2;
            double lim2 = 0;
            error = 10;
            double L = 2;
            double degree = 4;
            double speed = 3;
            double result2, result3;
            while(error > 0.000001)
            {
                step *= L;
                lim2 = lim1 + step;
                integral = NewtonKots(lim1,lim2,step);

                step = step / L;
                lim2 = lim1 + step;
                result2 = NewtonKots(lim1,lim2,step);

                step = step / L;
                lim2 = lim1 + step;
                result3 = NewtonKots(lim1,lim2,step);

                speed = -Math.Log(Math.Abs((result3 - result2) / (result2 - integral))) / Math.Log(L);
                System.Console.Out.Write(speed + "\n");

                error = Math.Abs((result2 - integral) / (Math.Pow(L, degree) - 1)); 
                integral += (result2 - integral) / (1 - Math.Pow(L, -degree));
            }
            WriteText("Составная ИКФ :",integral,Math.Abs(exactValue - integral),error);

            //* * * * * * * * * * * * * * * * Cоставная ИКФ c оптимальным шагом * * * * * * * * * *
            double hOpt = Math.Ceiling( (b - a) / (step * L * Math.Pow((0.00001 / error), 1.0 / degree)));
            hOpt = (b - a) / hOpt;
            lim1 = a;
            lim2 = a + hOpt;
            result2 = NewtonKots(lim1, lim2, hOpt);
            error = Math.Abs((result2 - integral) / (Math.Pow(L, degree) - 1)); 
            
            WriteText("Составная ИКФ с оптимальным шагом :",integral,Math.Abs(exactValue - result2),error);

            //* * * * * * * * * * * * * * * * * Квадратурная форма Гаусса * * * * * * * * * * * * * *

            lim1 = a;
            step = (b - a) / 10;
            error = 10;

            Console.WriteLine("скорость сходимости составной КФ Гаусса" + "\n");

            while(error > 1e-6)
            {
                
            }

        }

        static void WriteText(string s, double num)
            {
                System.Console.Out.Write(s + "\n");
                System.Console.Out.Write(num + "\n" + "\n");
            }

        static void WriteText(string s, double num, double num2, double num3)
        {
            System.Console.Out.Write("\n" + s + "\n");
            System.Console.Out.Write(num + "\n");
            System.Console.Out.Write("точная погрешность" + "\n");
            System.Console.Out.Write(num2 + "\n");
            System.Console.Out.Write("методическая погрешность" + "\n");
            System.Console.Out.Write(num3 + "\n" + "\n");
        }
    }
}
