using System;

namespace InfelicityTheory
{
    class Program
    {
        public static double Sqrt(double x)
        {
            const double eps = 1e-6/1.392;    //погрешность вычислений
            double xK=x;
            double xK1 = x;                   //число,корень которого нужно найти
            while(true)
            {
                xK1 = 0.5 * (xK + x / xK);    //следующий член
                if (Math.Abs(xK1 - xK) < eps)
                    break;
                xK = xK1;
            }
            return xK1;
        }
        public static double Ch(double x)
        {
            const double eps = 1e-6/1.392;   //погрешность вычислений
            double rowAmount = 0;             //сумма ряда
            double r = 1;                     //первый член ряда
            int rank = 0;                     //степень ряда
            while (Math.Abs(r) > eps)
            {
                rowAmount += r;
                rank=rank+1;
                r = r * x * x / ((2*rank) * (2*rank - 1));
            }
            return rowAmount+r;
        }

        public static double Sin(double x)
        {
            //while (x > Math.PI)
                //x -= Math.PI;
            const double eps = 1e-6/1.392;   //погрешность вычислений
            double rowAmount = 0;             //сумма ряда
            double r = x;                     //первый член ряда
            int rank = 1;                     //степень
            while(Math.Abs(r)>eps)
            {
                rowAmount += r;
                rank = rank + 2;              //наращиваем степень rank=1,3,5,...
                r = -r * x * x / (rank * (rank - 1));
            }
            return rowAmount+r;
        }
        static void Main(string[] args)
        {
            //Console.WriteLine(Program.Sin(0.2));
            //Console.WriteLine(Math.Sin(0.2));

            //Console.WriteLine(Program.Ch(0.2));
            //Console.WriteLine(Math.Cosh(0.2));
            //Console.WriteLine(Program.Sqrt(0.2));
            //Console.WriteLine(Math.Sqrt(0.2));
            double func;
            double defaultFunc;
            for(double i=0.2;i<0.3;i=i+0.01)
            {
                func = Program.Ch(Program.Sqrt(i * i + 0.3) / (1 + i)) * Program.Sin((1 + i) / (0.6 * i));
                defaultFunc = Math.Cosh(Math.Sqrt(i * i + 0.3) / (1 + i)) * Math.Sin((1 + i) / (0.6 * i));
                Console.WriteLine("{0}, {1} , {2}", func, defaultFunc,defaultFunc-func);
            }

        }
    }
}
