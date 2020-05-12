using System;

namespace Integral
{
    public class Functions
    {
        public double IntegralX0(double lim1,double lim2)
        {
            double sum;
            sum = 3 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 2 
                            - 3 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 2;
            return sum;
        }
        public double IntegralX1(double lim1,double lim2)
        {
            double sum;
            sum = 3.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 5 + 9.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 4 
                            - (3.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 5 + 9.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 4);
            return sum;
        } 
        public double IntegralX2(double lim1,double lim2)
        {
            double sum;
            sum = 3.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 8 + 9.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 5 + 27.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 8
                            - (3.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 8 + 9.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 5 + 27.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 8);
            return sum;
        }
        public double IntegralX3(double lim1,double lim2)
        {
            double sum;
            sum = 3.0 * Math.Pow((lim2 - 3.0 / 2),(11.0 / 3)) / 11 + 27.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 16 + 81.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 20 
                            + 81.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 16
                                - (3.0 * Math.Pow((lim1 - 3.0 / 2),(11.0 / 3)) / 11 + 27.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 16 + 81.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 20 
                                    + 81.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 16);
            return sum;
        }
        public double IntegralX4(double lim1,double lim2)
        {
            double sum;
            sum = 3.0 * Math.Pow((lim2 - 3.0 / 2),(14.0 / 3)) / 14 + 18.0 * Math.Pow((lim2 - 3.0 / 2),(11.0 / 3)) / 11 + 81.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 16 
                    + 81.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 10 + 243.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 32
                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(14.0 / 3)) / 14 + 18.0 * Math.Pow((lim1 - 3.0 / 2),(11.0 / 3)) / 11 + 81.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 16 
                    + 81.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 10 + 243.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 32);
            return sum;
        }
        public double IntegralX5(double lim1,double lim2)
        {
            double sum;
            sum = 3.0 * Math.Pow((lim2 - 3.0 / 2),(17.0 / 3)) / 17 + 45.0 * Math.Pow((lim2 - 3.0 / 2),(14.0 / 3)) / 28 + 135.0 * Math.Pow((lim2 - 3.0 / 2),(11.0 / 3)) / 22 
                    + 405.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 32 + 243.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 16 + 729.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 64
                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(17.0 / 3)) / 17 + 45.0 * Math.Pow((lim1 - 3.0 / 2),(14.0 / 3)) / 28 + 135.0 * Math.Pow((lim1 - 3.0 / 2),(11.0 / 3)) / 22 
                    + 405.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 32 + 243.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 16 + 729.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 64);
            return sum;
        }

        public double[] Kardano(double[] a, double[] x)
        {
             double current = a[0]; a[0] = a[2]; a[2] = current;
             double p = a[1] - a[0] * a[0] / 3.0; // p = b - a^2/3
             double q = a[2] + 2.0 * a[0] * a[0] * a[0] / 27.0 - a[0] * a[1] / 3.0; // q = c + 2a^3/27 - ab/3
             double determinant = q * q / 4.0 + p * p * p / 27.0;
             
             if (determinant < 0)
             {
                 double fi = 0;
                 if (q < 0) fi = Math.Atan(2.0 * Math.Sqrt(-determinant) / (-q));
                 if (q > 0) fi = Math.Atan(2.0 * Math.Sqrt(-determinant) / (-q) + Math.PI);
                 if (q == 0) fi = Math.PI / 2.0;

                 x[0] = 2.0 * Math.Sqrt(-p / 3.0) * Math.Cos(fi / 3.0) - a[0] / 3.0;
                 x[1] = 2.0 * Math.Sqrt(-p / 3.0) * Math.Cos(fi / 3.0 + 2.0 * Math.PI / 3.0) - a[0] / 3.0;
                 x[2] = 2.0 * Math.Sqrt(-p / 3.0) * Math.Cos(fi / 3.0 + 4.0 * Math.PI / 3.0) - a[0] / 3.0;
             }
             if (determinant > 0)
             {
                x[1] = 0;
                if ((-q) / 2.0 + Math.Pow(determinant, 1.0 / 2.0) < 0)
                    x[1] += -Math.Pow((q) / 2.0 - Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0);
                else x[1] += Math.Pow((-q) / 2.0 + Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0);
                if (-q / 2.0 - Math.Pow(determinant, 1.0 / 2.0) < 0)
                    x[1] += -Math.Pow(q / 2.0 + Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0) - a[0] / 3.0;
                else x[1] += Math.Pow(-q / 2.0 - Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0) - a[0] / 3.0;
             }
             if (determinant == 0){
                 x[0] = 2 * Math.Pow(-q / 2.0, 1.0 / 3.0) - a[0] / 3.0;
                 x[1] =  -Math.Pow(-q / 2.0, 1.0 / 3.0) - a[0] / 3.0;
             }
            return x;
        }
    }
}