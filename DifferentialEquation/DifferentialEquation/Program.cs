﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.Remoting.Messaging;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace DifferentialEquation
{
    class Program
    {
        delegate Complex Function(Complex x, Complex[] y);

        // Двухэтапный ЯМРК второго порядка, зависящий от константы с2
        static Complex[,] RKMethod(Function[] derivatives, double steps, double c2, Complex[] y0, double xStart, double xEnd)
        {
            var h = (xEnd - xStart) / steps;
            var mas = new Complex[4, (int)steps + 1];

            var x = xStart;
            var y = y0;

            for (int i = 1; i <= steps; i++)
            {

                y = IterationMethod(derivatives, h, c2, y, x);
                for (int j = 0; j < 4; j++)
                    mas[j, i] = y[j];
                x += h;

            }

            return mas;
        }


        // Метод с автоматическим выбором шага
        static Solution AutomaticMethod(Function[] derivatives, double h1, double c2, Complex[] y0, double xStart, double xEnd, double rtol, double atol)
        {
            var h = h1;
            var x = xStart;
            var maxStep = 0.01;
            var y = y0;

            int counter = 0;
            var rn = new Complex[4];
            var rn2 = new Complex[4];
            do
            {
                counter++;
                var yf = IterationMethod(derivatives, h, c2, y, x);
                var y1 = IterationMethod(derivatives, h / 2, c2, y, x);
                var y2 = IterationMethod(derivatives, h / 2, c2, y1, x + h / 2);

                for (int i = 0; i < 4; i++)
                {
                    rn[i] = (y2[i].Real - yf[i].Real) / (1.0 - Math.Pow(2, -2));
                    rn2[i] = (y2[i].Real - yf[i].Real) / (Math.Pow(2, 2) - 1.0);
                }

                var rnorm = norma(rn);
                var yfnorm = norma(yf);

                if (rnorm > ((yfnorm * rtol + atol)) * 4)
                {
                    h /= 2;
                    //y = y0;
                }
                else if ((yfnorm * rtol + atol) < rnorm)
                {
                    x += h;
                    h /= 2;
                    for (int i = 0; i < 4; i++)
                        y[i] = y2[i] + rn2[i];
                }
                else if ((yfnorm * rtol + atol) / 8.0 <= rnorm)
                {
                    x += h;
                    for (int i = 0; i < 4; i++)
                        y[i] = yf[i] + rn[i];
                }
                else
                {
                    x += h;
                    h = 2 * h;
                    for (int i = 0; i < 4; i++)
                        y[i] = yf[i] + rn[i];
                }

                if (h > maxStep) h = maxStep;

                if ((x < xEnd) && ((x + h) > xEnd))
                {
                    h = xEnd - x;
                }

            } while (x < xEnd);

            Solution sol = new Solution();
            sol.count = counter;
            sol.answer = y;
            return sol;

            //System.Console.Out.WriteLine("Число шагов автоматического метода " + counter);
            //return y;
        }

        // Евклидова норма
        static double norma(Complex[] v)
        {
            double n = 0;
            for (int i = 0; i < v.Length; i++)
                n += Math.Pow(v[i].Real, 2);
            return Math.Sqrt(n);
        }

        // Итерация двухэтапного метода Рунге-Кутты второго порядка, зависящая от константы с
        static Complex[] IterationMethod(Function[] derivatives, double h, double c, Complex[] y0, double x0)
        {
            double[] b = new double[2];

            Complex[] V1 = new Complex[4];
            Complex[] V2 = new Complex[4];
            Complex[] Vt = new Complex[4];

            // conditions of order
            b[1] = 1 / (2 * c);
            b[0] = 1 - b[1];

            for (int i = 0; i < 4; i++)
                V1[i] = derivatives[i](x0, y0);

            for (int i = 0; i < 4; i++)
                Vt[i] = y0[i] + c * h * V1[i];

            for (int i = 0; i < 4; i++)
                V2[i] = derivatives[i](x0 + c * h, Vt);

            var F = new Complex[4];

            for (int p = 0; p < 4; p++)
            {
                F[p] = y0[p] + h * (b[0] * V1[p] + b[1] * V2[p]);
            }
            return F;
        }

        static void output(double mas1, double mas2, double mas3, double mas4)
        {
            Console.WriteLine("x=5: y0 = " + mas1);
            Console.WriteLine("x=5: y1 = " + mas2);
            Console.WriteLine("x=5: y2 = " + mas3);
            Console.WriteLine("x=5: y3 = " + mas4 + "\n");
        }

        public class Solution
        {
            public Complex[] answer;
            public double error;
            public int count;
        }

        
        static Solution AutomiticSolution(double rtol, Function[] Derivatives, Complex[] y0, double acc, double acc1, double acc2, double acc3, double c2)
        {
            double atol = 1e-12;
            var norm = Math.Pow((Derivatives[0](0, y0) * Derivatives[0](0, y0)).Real + (Derivatives[1](0, y0) * Derivatives[1](0, y0)).Real
                            + (Derivatives[2](0, y0) * Derivatives[2](0, y0)).Real + (Derivatives[3](0, y0) * Derivatives[3](0, y0)).Real, 1 / 2.0);

            var delta = Math.Pow(1 / 5.0, 3) + Math.Pow(norm, 3);
            var h1 = Math.Pow((rtol * norm + atol) / delta, 1.0 / 3);
            
            Solution sol = AutomaticMethod(Derivatives, h1, c2, y0, 0, 5, rtol, atol);
            var autoError = Math.Sqrt(Math.Pow(sol.answer[0].Real - acc, 2) + Math.Pow(sol.answer[1].Real - acc1, 2)
                            + Math.Pow(sol.answer[2].Real - acc2, 2) + Math.Pow(sol.answer[3].Real - acc3, 2));

            sol.error = autoError;

            return sol;
        }

        static void Main(string[] args)
        {
            double A = -3.0;
            double B = 2.0;
            double C = 1.0; 

            // Initial derivatives
            var Derivatives = new Function[4];
            Derivatives[0] = (x, y) => (2 * x * Complex.Pow(y[1], 1.0 / B) * y[3]);
            Derivatives[1] = (x, y) => (2 * B * x * Complex.Exp((B / C)*(y[2] - A)) * y[3]);
            Derivatives[2] = (x, y) => (2 * C * x * y[3]);
            Derivatives[3] = (x, y) => (-2 * x * Complex.Log(y[0]));

            // Initial conditions of y1,y2,y3,y4
            var y0 = new Complex[] {1,1,A,1}; // 

            var acc = Math.Exp(Math.Sin(5 * 5));
            var acc1 = Math.Exp(B * Math.Sin(5 * 5));
            var acc2 = (C * Math.Sin(5 * 5) + A);
            var acc3 = Math.Cos(5 * 5);

            Console.WriteLine("Точные значения");
            output(acc, acc1, acc2, acc3);

            int steps = 1300; //Для удобства вывода графиков
            
            var mas = RKMethod(Derivatives, 1300, 0.1, y0, 0, 5);
            Console.WriteLine("Двухэтапный ЯМРК второго порядка, с2 = 0.1");
            Console.WriteLine("1300 шагов, x=5: y0 = " + mas[0, steps].Real);

            var grid2 = RKMethod(Derivatives, steps, 1, y0, 0, 5);
            Console.WriteLine("Метод Хойна, с2 = 1");
            Console.WriteLine("1300 шагов, x=5: y0 = " + grid2[0, steps].Real);
            Console.WriteLine("");

            double[] mist = new double[6];
            for (int i = 1; i <= 6; i++)
            {
                var masl = RKMethod(Derivatives, 5*Math.Pow(2, i), 0.1, y0, 0, 5);
                mist[i-1] = Math.Sqrt(Math.Pow(masl[0, (int)(5 * Math.Pow(2, i))].Real - Math.Exp(Math.Sin(5 * 5)), 2)
                    + Math.Pow(masl[1, (int)(5 * Math.Pow(2, i))].Real - Math.Exp(B*Math.Sin(5 * 5)), 2)
                    + Math.Pow(masl[2, (int)(5 * Math.Pow(2, i))].Real - (C * Math.Sin(5 * 5) + A), 2)
                    + Math.Pow(masl[3, (int)(5 * Math.Pow(2, i))].Real - Math.Cos(5 * 5), 2)
                    );
            }
            DrawGraph.Draw(0, 0, mas, grid2, mist);

            // Таблица "Порядок сходимости методов"
            Console.WriteLine("Двухэтапный ЯМРК второго порядка, с2 = 0.1, х=5, результаты для у0");

            double h = 0.01;
            var er_pred = 0.08;
            Console.WriteLine("{0,12} | {1, 10}|  {2,20} | {3,14}", "Число шагов", "h", "Error", "log2(err/er)");
            for (int i = 0; i <= 6; i++)
            {
                if (i != 0) er_pred = Math.Abs(mas[0, steps].Real - mas[0, steps - 1].Real);

                h /= 2;
                steps = (int) (5.0 / h);
                mas = RKMethod(Derivatives, steps, 0.1, y0, 0, 5);
                Console.WriteLine("{0,12} | {1, 10}|  {2,20} | {3,14}", steps, h, acc - mas[0, steps].Real,
                    Math.Log(Math.Pow(er_pred / Math.Abs(mas[0, steps].Real - mas[0, steps - 1].Real), 2.0), 2.0));
            } 
            Console.WriteLine("");

            // Таблица для метода Хойна
            Console.WriteLine("Метод хойна, с2 = 1, х=5, результаты для у0");

            h = 0.01;
            er_pred = 0.08;
            Console.WriteLine("{0,12} | {1, 10}|  {2,20} | {3,14}", "Число шагов", "h", "Error", "log2(err/er)");
            for (int i = 0; i <= 6; i++)
            {
                if (i != 0) er_pred = Math.Abs(mas[0, steps].Real - mas[0, steps-1].Real);
                
                h /= 2;
                steps = (int) (5.0 / h);
                mas = RKMethod(Derivatives, steps, 1, y0, 0, 5);

                Console.WriteLine("{0,12} | {1, 10}|  {2,20} | {3,14}", steps, h, acc - mas[0, steps].Real,
                    Math.Log(Math.Pow(er_pred / Math.Abs(mas[0, steps].Real - mas[0, steps-1].Real), 2.0), 2.0));

            }
            Console.WriteLine("");

            // метод с выбором оптимального шага
            var mash = RKMethod(Derivatives, 5.0 * 50, 0.1, y0, 0, 5);
            var mash2 = RKMethod(Derivatives, 5.0 * 100, 0.1, y0, 0, 5);

            var Rn = new Complex[4];

            for (int i = 0; i < 4; i++)
                Rn[i] = (mash2[i, (int) (5.0*100)].Real - mash[i, (int) (5.0*50)].Real)/(1 - Math.Pow(2, -2));

            double tol = 1e-6;
            var hopt = (1.0/50)*Math.Pow(tol/(Math.Abs( norma(Rn) )), 1/2.0);
            var masht = RKMethod(Derivatives, (int)(Math.Ceiling(5 / hopt)), 0.1, y0, 0, 5);

            var error = Math.Sqrt(Math.Pow(masht[0, (int)(Math.Ceiling(5 / hopt))].Real - acc, 2) + Math.Pow(masht[1, (int)(Math.Ceiling(5 / hopt))].Real - acc1, 2)
                            + Math.Pow(masht[2, (int)(Math.Ceiling(5 / hopt))].Real - acc2, 2) + Math.Pow(masht[3, (int)(Math.Ceiling(5 / hopt))].Real - acc3, 2));

            Console.WriteLine("Оптимальный шаг:" + hopt);
            Console.WriteLine("Число шагов:" + Math.Ceiling(5 / hopt));
            Console.WriteLine("Ошибка x=5: " + error);
            //output(masht[0, (int)(Math.Ceiling(5 / hopt))].Real, masht[1, (int)(Math.Ceiling(5 / hopt))].Real,
            //masht[2, (int)(Math.Ceiling(5 / hopt))].Real, masht[3, (int)(Math.Ceiling(5 / hopt))].Real);
            Console.WriteLine("");


            Solution sol = new Solution();

            // Автоматический метод из условия, с2 = 0.1
            Console.WriteLine("ЯМРК с2 = 0.1 c автоматическим выбором шага");
            double rtol = 1e-6;
            double c2 = 0.1;
            sol = AutomiticSolution(rtol, Derivatives, y0, acc, acc1, acc2, acc3, c2);

            Console.WriteLine("Ошибка: " + sol.error);
            Console.WriteLine("Количество вычислений: " + sol.count);
            //Console.WriteLine("Ответ х=5: у0 = " + sol.answer[0].Real + "\n");
            Console.WriteLine("");


            // Автоматический метод Хойна, с2 = 1
            Console.WriteLine("Метода Хойна (с2 = 1) с автоматическим выбором шага ");
            rtol = 1e-6;
            c2 = 1;
            sol = AutomiticSolution(rtol, Derivatives, y0, acc, acc1, acc2, acc3, c2);

            Console.WriteLine("Ошибка: " + sol.error);
            Console.WriteLine("Количество вычислений: " + sol.count);
            //Console.WriteLine("Ответ х=5: у0 = " + sol.answer[0].Real + "\n");

            // Таблица Затраты-Точность
            Console.WriteLine("Таблица затраты - точность");
            Console.WriteLine("{0,12} | {1, 20}|  {2,10} ", "tol", "Error", "Count");
            rtol = 1e-6;
            c2 = 0.1;
            for (int i = 0; i < 6; i++)
            {
                sol = AutomiticSolution(rtol, Derivatives, y0, acc, acc1, acc2, acc3, c2);
                Console.WriteLine("{0,12} | {1, 20}|  {2,10} ", rtol, sol.error, sol.count);
                rtol /= 10;
            }

            Console.Read();

        }

       
    }
}
