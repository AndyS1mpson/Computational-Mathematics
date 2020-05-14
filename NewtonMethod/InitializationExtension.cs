using System;

namespace NewtonMethod
{
    public static class InitializationExtension
    {
        // Initialize function 
        public static void InitF(this double[] F,double[] x)
        {
            F[0] = Math.Sin(x[0] + 1) - x[1] - 1.2;
            F[1] = 2 * x[0] + Math.Cos(x[1]) - 2;
        }
        public static double[,] InitJ(double[] x)
        {
            double[,] matr = new double[2,2];
            matr[0,0] = Math.Cos(x[0] + 1);
            matr[0,1] = -1;
            matr[1,0] = 2;
            matr[1,1] = -Math.Sin(x[1]);
            return matr;
        }
        public static void InitX(this double[] X)
        {
            // x[0] = 1;
            // x[1] = 2;
            X[0] = 0.5; X[1] = 0.5; X[2] = 1.5; X[3] = -1; X[4] = -0.2; 
            X[5] = 1.5; X[6] = 0.5; X[7] = -0.5; X[8] = 1.5; X[9] = -1.5;
        }
    
        public static void InitializationF(this double[] F,double[] X) {        
            F[0] = Math.Cos(X[0] * X[1]) - Math.Exp(-3.0 * X[2]) + X[3] * X[4] * X[4] - X[5] - Math.Sinh(2.0 * X[7]) * X[8] + 2.0 * X[9] + 2.0004339741653854440;
            F[1] = Math.Sin(X[0] * X[1]) + X[2] * X[8] * X[6] - Math.Exp(-X[9] + X[5]) + 3.0 * X[4] * X[4] - X[5] * (X[7] + 1.0) + 10.886272036407019994;
            F[2] = X[0] - X[1] + X[2] - X[3] + X[4] - X[5] + X[6] - X[7] + X[8] - X[9] - 3.1361904761904761904;
            F[3] = 2.0 * Math.Cos(-X[8] + X[3]) + X[4] / (X[2] + X[0]) - Math.Sin(X[1] * X[1]) + Math.Cos(X[6] * X[9]) * Math.Cos(X[6] * X[9]) - X[7] - 0.1707472705022304757;
            F[4] = Math.Sin(X[4]) + 2.0 * X[7] * (X[2] + X[0]) - Math.Exp(-X[6] * (-X[9] + X[5])) + 2.0 * Math.Cos(X[1]) - 1.00 / (X[3] - X[8]) - 0.3685896273101277862;
            F[5] = Math.Exp(X[0] - X[3] - X[8]) + X[4] * X[4] / X[7] + Math.Cos(3.0 * X[9] * X[1]) / 2.0 - X[5] * X[2] + 2.0491086016771875115;
            F[6] = X[1] * X[1] * X[1] * X[6] - Math.Sin(X[9] / X[4] + X[7]) + (X[0] - X[5]) * Math.Cos(X[3]) + X[2] - 0.7380430076202798014;
            F[7] = X[4] * (X[0] - 2.0 * X[5]) * (X[0] - 2.0 * X[5]) - 2.0 * Math.Sin(-X[8] + X[2]) + 1.5 * X[3] - Math.Exp(X[1] * X[6] + X[9]) + 3.5668321989693809040;
            F[8] = 7.0 / X[5] + Math.Exp(X[4] + X[3]) - 2.0 * X[1] * X[7] * X[9] * X[6] + 3.0 * X[8] - 3.0 * X[0] - 8.4394734508383257499;
            F[9] = X[9] * X[0] + X[8] * X[1] - X[7] * X[2] + Math.Sin(X[3] + X[4] + X[5]) * X[6] - 0.78238095238095238096;
        }

        public static double[,] InitializationJ(double[] X)
        {
            double[,] matrix = new double[10,10];
            matrix[0,0] = -Math.Sin(X[0] * X[1]) * X[1];
            matrix[0,1] = -Math.Sin(X[0] * X[1]) * X[0];
            matrix[0,2] = 3 * Math.Exp(-(3 * X[2]));
            matrix[0,3] = X[4] * X[4];
            matrix[0,4] = 2 * X[3] * X[4];
            matrix[0,5] = -1;
            matrix[0,6] = 0;
            matrix[0,7] = -2 * Math.Cosh((2 * X[7])) * X[8];
            matrix[0,8] = -Math.Sinh((2 * X[7]));
            matrix[0,9] = 2;
            matrix[1,0] = Math.Cos(X[0] * X[1]) * X[1];
            matrix[1,1] = Math.Cos(X[0] * X[1]) * X[0];
            matrix[1,2] = X[8] * X[6];
            matrix[1,3] = 0;
            matrix[1,4] = 6 * X[4];
            matrix[1,5] = -Math.Exp(-X[9] + X[5]) - X[7] - 0.1e1;
            matrix[1,6] = X[2] * X[8];
            matrix[1,7] = -X[5];
            matrix[1,8] = X[2] * X[6];
            matrix[1,9] = Math.Exp(-X[9] + X[5]);
            matrix[2,0] = 1;
            matrix[2,1] = -1;
            matrix[2,2] = 1;
            matrix[2,3] = -1;
            matrix[2,4] = 1;
            matrix[2,5] = -1;
            matrix[2,6] = 1;
            matrix[2,7] = -1;
            matrix[2,8] = 1;
            matrix[2,9] = -1;
            matrix[3,0] = -X[4] * Math.Pow(X[2] + X[0], -2);
            matrix[3,1] = -2 * Math.Cos(X[1] * X[1]) * X[1];
            matrix[3,2] = -X[4] * Math.Pow(X[2] + X[0], -2);
            matrix[3,3] = -2 * Math.Sin(-X[8] + X[3]);
            matrix[3,4] = 1 / (X[2] + X[0]);
            matrix[3,5] = 0;
            matrix[3,6] = -2 * Math.Cos(X[6] * X[9]) * Math.Sin(X[6] * X[9]) * X[9];
            matrix[3,7] = -1;
            matrix[3,8] = 2 * Math.Sin(-X[8] + X[3]);
            matrix[3,9] = -2 * Math.Cos(X[6] * X[9]) * Math.Sin(X[6] * X[9]) * X[6];
            matrix[4,0] = 2 * X[7];
            matrix[4,1] = -2 * Math.Sin(X[1]);
            matrix[4,2] = 2 * X[7];
            matrix[4,3] = Math.Pow(-X[8] + X[3], -2);
            matrix[4,4] = Math.Cos(X[4]);
            matrix[4,5] = X[6] * Math.Exp(-X[6] * (-X[9] + X[5]));
            matrix[4,6] = -(X[9] - X[5]) * Math.Exp(-X[6] * (-X[9] + X[5]));
            matrix[4,7] = (2 * X[2]) + 2 * X[0];
            matrix[4,8] = -Math.Pow(-X[8] + X[3], -2);
            matrix[4,9] = -X[6] * Math.Exp(-X[6] * (-X[9] + X[5]));
            matrix[5,0] = Math.Exp(X[0] - X[3] - X[8]);
            matrix[5,1] = -3.0 / 2.0 * Math.Sin(3 * X[9] * X[1]) * X[9];
            matrix[5,2] = -X[5];
            matrix[5,3] = -Math.Exp(X[0] - X[3] - X[8]);
            matrix[5,4] = 2 * X[4] / X[7];
            matrix[5,5] = -X[2];
            matrix[5,6] = 0;
            matrix[5,7] = -X[4] * X[4] *  Math.Pow(X[7], (-2));
            matrix[5,8] = -Math.Exp(X[0] - X[3] - X[8]);
            matrix[5,9] = -3.0 / 2.0 * Math.Sin(3 * X[9] * X[1]) * X[1];
            matrix[6,0] = Math.Cos(X[3]);
            matrix[6,1] = 3 * X[1] * X[1] * X[6];
            matrix[6,2] = 1;
            matrix[6,3] = -(X[0] - X[5]) * Math.Sin(X[3]);
            matrix[6,4] = Math.Cos(X[9] / X[4] + X[7]) * X[9] * Math.Pow(X[4], (-2));
            matrix[6,5] = -Math.Cos(X[3]);
            matrix[6,6] = Math.Pow(X[1], 3);
            matrix[6,7] = -Math.Cos(X[9] / X[4] + X[7]);
            matrix[6,8] = 0;
            matrix[6,9] = -Math.Cos(X[9] / X[4] + X[7]) / X[4];
            matrix[7,0] = 2 *  X[4] * (X[0] - 2 * X[5]);
            matrix[7,1] = -X[6] * Math.Exp(X[1] * X[6] + X[9]);
            matrix[7,2] = -2 * Math.Cos(-X[8] + X[2]);
            matrix[7,3] = 0.15e1;
            matrix[7,4] = Math.Pow(X[0] - 2 * X[5], 2);
            matrix[7,5] = -4 *  X[4] * (X[0] - 2 * X[5]);
            matrix[7,6] = -X[1] * Math.Exp(X[1] * X[6] + X[9]);
            matrix[7,7] = 0;
            matrix[7,8] = 2 * Math.Cos(-X[8] + X[2]);
            matrix[7,9] = -Math.Exp(X[1] * X[6] + X[9]);
            matrix[8,0] = -3;
            matrix[8,1] = -2 *  X[7] * X[9] * X[6];
            matrix[8,2] = 0;
            matrix[8,3] = Math.Exp((X[4] + X[3]));
            matrix[8,4] = Math.Exp((X[4] + X[3]));
            matrix[8,5] = -0.7e1 * Math.Pow(X[5], -2);
            matrix[8,6] = -2 * X[1] *  X[7] * X[9];
            matrix[8,7] = -2 * X[1] * X[9] * X[6];
            matrix[8,8] = 3;
            matrix[8,9] = -2 * X[1] *  X[7] * X[6];
            matrix[9,0] = X[9];
            matrix[9,1] = X[8];
            matrix[9,2] = -X[7];
            matrix[9,3] = Math.Cos(X[3] + X[4] + X[5]) * X[6];
            matrix[9,4] = Math.Cos(X[3] + X[4] + X[5]) * X[6];
            matrix[9,5] = Math.Cos(X[3] + X[4] + X[5]) * X[6];
            matrix[9,6] = Math.Sin(X[3] + X[4] + X[5]);
            matrix[9,7] = -X[2];
            matrix[9,8] = X[1];
            matrix[9,9] = X[0];

            return matrix;
        }
    }
}