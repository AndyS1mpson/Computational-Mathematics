using System;

namespace Integral
{
    class Program
    {
        // my function 
        static double f(double x)
            => 2 * Math.Cos(2.5 * x) * Math.Exp(x * 1.0 / 3) + 4 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + x;
        static double NewtonKotsKF(double lim1, double lim2, double step)
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
        
        static double GaussKF(double lim1,double lim2,double step)
        {
            double result = 0;
            while(lim2 <= b)
            {
                double[] moments = new double[6];
                moments[0] = functions.IntegralX0(lim1,lim2);
                moments[1] = functions.IntegralX1(lim1,lim2);
                moments[2] = functions.IntegralX2(lim1,lim2);
                moments[3] = functions.IntegralX3(lim1,lim2);
                moments[4] = functions.IntegralX4(lim1,lim2);
                moments[5] = functions.IntegralX5(lim1,lim2);


                // add 3 new nodes
                double[] nodes = new double[3];
                nodes[0] = lim1;
                nodes[1] = lim1 + (lim2-lim1)/2;
                nodes[2] = lim2;
                // declare the coeff of system 
                double[,] matrixArray = new double[3,3];
                matrixArray[0,0] = moments[0];
                matrixArray[0,1] = moments[1];
                matrixArray[0,2] = moments[2];
                matrixArray[1,0] = moments[1];
                matrixArray[1,1] = moments[2];
                matrixArray[1,2] = moments[3];
                matrixArray[2,0] = moments[2];
                matrixArray[2,1] = moments[3];
                matrixArray[2,2] = moments[4];
                
                // declare matrix
                SquareMatrix matrix = new SquareMatrix(matrixArray,3);
                matrix.LUDecomposition();

                double[] y = new double[3];
                y[0] = -moments[3];
                y[1] = -moments[4];
                y[2] = -moments[5];

                // find the solution
                var solution = matrix.SolutionSystem(y);
                double[] solv = new double[solution.GetLength(1)];
                for(int i = 0;i < solution.GetLength(1);i++)
                    solv[i] = solution[0,i];

                var p = functions.Kardano(solv,nodes);

                Array.Sort(p);

                matrix._matrix[0,0] = 1;
                matrix._matrix[0,1] = 1;
                matrix._matrix[0,2] = 1;
                matrix._matrix[1,0] = p[0];
                matrix._matrix[1,1] = p[1];
                matrix._matrix[1,2] = p[2];
                matrix._matrix[2,0] = p[0] * p[0];
                matrix._matrix[2,1] = p[1] * p[1];
                matrix._matrix[2,2] = p[2] * p[2];

                matrix.LUDecomposition();
                
                var s = matrix.SolutionSystem(moments);
                double[] A = new double[s.GetLength(1)];
                for(int i = 0;i < A.Length;i++)
                    A[i] = s[0,i];

                for (int i = 0; i < 3; i++)
                {
                    result += A[i] * f(nodes[i]);
                }
                lim1 = lim2;
                lim2 += step;

            }
            return result;
        }
        static double a = 1.5;
        static    double b = 3.3;
        static    double alpha = 1/3;
        static    double beta = 0;
        static    double exactValue = 7.077031437995793610263911711602477164432;
        static    Functions functions = new Functions();
        static void Main(string[] args)
        {

            // * * * * * * * * * * * * * Вариант Ньютона - Котса * * * * * * * * * * * *
            Console.WriteLine("Вариант Ньютона-Котса : " + "\n");

            
            double integral = NewtonKotsKF(a,b,1);
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
            Console.WriteLine("Скорость сходимости составной ИКФ :");
            while(error > 0.000001)
            {
                step *= L;
                lim2 = lim1 + step;
                integral = NewtonKotsKF(lim1,lim2,step);

                step = step / L;
                lim2 = lim1 + step;
                result2 = NewtonKotsKF(lim1,lim2,step);

                step = step / L;
                lim2 = lim1 + step;
                result3 = NewtonKotsKF(lim1,lim2,step);

                speed = -Math.Log(Math.Abs((result3 - result2) / (result2 - integral))) / Math.Log(L);
                Console.WriteLine(speed);

                error = Math.Abs((result2 - integral) / (Math.Pow(L, degree) - 1)); 
                integral += (result2 - integral) / (1 - Math.Pow(L, -degree));
            }
            Console.WriteLine("\n");
            WriteText("Составная ИКФ :",integral,Math.Abs(exactValue - integral),error);

            //* * * * * * * * * * * * * * * * Cоставная ИКФ c оптимальным шагом * * * * * * * * * *
            double hOpt = Math.Ceiling( (b - a) / (step * L * Math.Pow((0.00001 / error), 1.0 / degree)));
            hOpt = (b - a) / hOpt;
            lim1 = a;
            lim2 = a + hOpt;
            result2 = NewtonKotsKF(lim1, lim2, hOpt);
            error = Math.Abs((result2 - integral) / (Math.Pow(L, degree) - 1)); 
            
            WriteText("Составная ИКФ с оптимальным шагом :",integral,Math.Abs(exactValue - result2),error);

            //* * * * * * * * * * * * * * * * * Квадратурная форма Гаусса * * * * * * * * * * * * * *

            lim1 = a;
            step = (b - a) / 10;
            error = 10;

            Console.WriteLine("Скорость сходимости составной КФ Гаусса" + "\n");

            while(error > 1e-6)
            {
                 step *= L;
                lim2 = a + step;
                integral = GaussKF(lim1, lim2, step);

                step = step / L;
                lim2 = a + step;
                result2 = GaussKF(lim1, lim2, step);

                step = step / L ;
                lim2 = a + step;
                result3 = GaussKF(lim1, lim2, step);

                speed = -Math.Log(Math.Abs((result3 - result2) / (result2 - integral))) / Math.Log(L);
                System.Console.Out.Write(speed + "\n");
                error = Math.Abs((result2 - integral) / (Math.Pow(L, degree) - 1));
                integral += (result2 - integral) / (1 - Math.Pow(L, -degree));
            }

            WriteText("Составная ИКФ по Гауссу", integral, Math.Abs(exactValue - integral), error);

             //* * * * * * * * * * * * * * * Cоставная форма Гаусса c оптимальным шагом * * * * * * * * * *
            hOpt = Math.Ceiling((b - a) / (step * L * Math.Pow((0.00001 / error), 1.0 / degree)));
            hOpt = (b - a) / hOpt;
            lim1 = a;
            lim2 = a + hOpt;
            result2 = GaussKF(lim1, lim2, hOpt);
            error = Math.Abs((result2 - integral) / (Math.Pow(L, degree) - 1));

            WriteText("составная ИКФ Гаусса с оптимальным шагом", result2, Math.Abs(exactValue - result2), error);
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
