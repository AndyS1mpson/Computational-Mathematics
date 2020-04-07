using System;
using System.Diagnostics;

namespace NewtonMethod
{
    public class NewtonMethod
    {
        private double[] F;                     // system of functions
        private double[] x,xk1;                 // approximate roots
        private double[] currentX;              
        private SquareMatrix J, I;              // matrix of Jacoby 
        private int size;                       // size of matrices
        private double eps;
        private double norm;                    // norm of vector
        private int iter;                       // amount of iterations
        public Stopwatch sw;                    // check time of work
        public NewtonMethod()
        { 
            //size = 2;
            size = 10;
            eps = 1e-3;
            x = new double[size];
            xk1 = new double[size];
            currentX = new double[size];
            F = new double[size];
            iter = 0;
            norm = 0;
            sw = new Stopwatch();
        }

        // Newton's method through the search of inverse matrix
        public double[] Method()
        {
            Console.Out.Write("Метод Ньютона с поиском обратной матрицы :");
            sw.Start();


            x.InitX();
            iter = 0;
            while(true)
            {
                iter++;
                // initialize function F and Jacoby's matrix
                F.InitializationF(x);
                norm = 0;
                J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);
                // make LU decomposition
                J.LUDecomposition();
                
                currentX = SquareMatrix.multiplyVec(J.reverseMatrix,F,size);
                
                for(int i = 0;i < size;i++)
                {
                    xk1[i] = x[i] - currentX[i];

                    norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);

                    x[i] = xk1[i];
                }

                norm = Math.Sqrt(norm);
                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");
            return xk1;
        }

        // Modified Newton's method with LU decomposition
        public double[] ModifiedMethod()
        {
           Console.Out.Write("Модифицированный метод Ньютона :");
           sw.Start();
           x.InitX();
           iter = 0;
           //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
           J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);

           J.LUDecomposition();
           J.Reverse();

            while(true)
            {
               iter++;
               norm = 0;
               //F.InitF(x);
               F.InitializationF(x);
               currentX = J.SolutionSystem(F);

                for(int i = 0;i < size;i++)
                {
                   xk1[i] = x[i] - currentX[i];

                   norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);

                    x[i] = xk1[i];
                }
                
                norm = Math.Sqrt(norm);
                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");

            return xk1;
         }   
    
        // Method with transition to modified 
        public double[] TransitionToModMeth()
        {
            Console.Out.Write("Переход от обычного метода Ньютона к модифицированному :");
            int step = 7;

            x.InitX();

            iter = 0;
            sw.Start();
            while(true)
            {
                iter++;
                norm = 0;
                //F.InitF(x);
                //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
                J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);
                F.InitializationF(x);
                J.LUDecomposition();
                if(step > 0)
                {
                    J.Reverse();
                    currentX = SquareMatrix.multiplyVec(J.reverseMatrix,F,size);
                    step--;
                }
                else
                {
                    currentX = J.SolutionSystem(F);
                }
                for(int i = 0;i < size;i++)
                {
                    xk1[i] = x[i] - currentX[i];
                    norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);

                    x[i] = xk1[i];
                }
                 norm = Math.Sqrt(norm);
                if(norm < eps)
                    break;

            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");
            return xk1;
        }

        // method with recounting reverse matrix every k iterations
        public double[] MethodWithRewriteRevMat()
        {
            int step = 5;
            Console.Out.Write("Поиск обратной матрицы каждые k = " + step + " итераций :");
            
            
            x.InitX();
            //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
            J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);

            J.LUDecomposition();
            iter = 0;
            
            sw.Start();
            while(true)
            {
                iter++;
                norm = 0;
                //F.InitF(x);
                F.InitializationF(x);
                if(iter % step == 1)
                {
                    //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
                J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);

                    J.LUDecomposition();
                    J.Reverse();
                }
                currentX = J.SolutionSystem(F);

                for(int i = 0;i < size;i++)
                {
                    xk1[i] = x[i] - currentX[i];

                    norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);

                    x[i] = xk1[i];
                }

                norm = Math.Sqrt(norm);
                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");
            return xk1;
        }
    }
}