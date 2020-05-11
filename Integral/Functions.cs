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
    }
}