using System;

namespace Integral
{
    class Program
    {
        static void Main(string[] args)
        {
            double a = 1.5;
            double b = 3.3;
            double alpha = 1/3;
            double beta = 0;

            double[] nodes;
            double[] moments;
            // add 3 new nodes
            nodes = new double[3];
            nodes[0] = a + (b - a) / 4;
            nodes[1] = a + 2 * (b - a) / 4;
            nodes[2] = a + 3 * (b - a) / 4;

            // find out moments
            moments = new double[3];
            moments[0] = 2.21959;
            moments[1] = 4.92749;
            moments[2] = 11.5863;

            double[] matrix = new double[9];т
        }
    }
}
