using System;
using System.Diagnostics;

namespace NewtonMethod
{
    public class NewtonMethod
    {
        private double[,] F;                     // system of functions
        private double[] F1;
        private double[] x,xk1;                 // approximate roots
        private double[] currentX;              
        private SquareMatrix J;              // matrix of Jacoby 
        private int size;                       // size of matrices
        private double eps;
        private int numOfOperations;
        private double norm;                    // norm of vector
        private int iter;                       // amount of iterations
        public Stopwatch sw;                    // check time of work
        public NewtonMethod()
        { 
            //size = 2;
            size = 10;
            F = new double[1,size];
            eps = 1e-3;
            x = new double[size];
            xk1 = new double[size];
            currentX = new double[size];
            F1 = new double[size];
            iter = 0;
            norm = 0;
            sw = new Stopwatch();
        }

        // Newton's method through the search of inverse matrix
        public double[] Method()
        {
            Console.Out.Write("Метод Ньютона с поиском обратной матрицы :" + "\n");
            sw.Start();


            x.InitX();
            iter = 0;
            numOfOperations = 0;
            while(true)
            {
                iter++;

                // initialize function F and Jacoby's matrix
                F1.InitializationF(x);
                for(int i = 0;i < size;i++)
                    F[0,i] = F1[i];
                norm = 0;
                J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);
                // make LU decomposition
                J.LUDecomposition();
                numOfOperations += J.numOfOperations;

                currentX = SquareMatrix.multiplyVec(J.reverseMatrix,F1,size);
                
                for(int i = 0;i < size;i++)
                {
                    xk1[i] = x[i] - currentX[i];

                    //norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);
                
                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);


                    x[i] = xk1[i];

                    numOfOperations += 5;
                }

                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Число арифметических операций :" + numOfOperations + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");
            return xk1;
        }

        // Modified Newton's method with LU decomposition
        public double[] ModifiedMethod()
        {
           Console.Out.Write("Модифицированный метод Ньютона :" + "\n");
           sw.Start();
           x.InitX();
           iter = 0;
           numOfOperations = 0;
           //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
           J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);

           J.LUDecomposition();
           numOfOperations += J.numOfOperations;
           //J.Reverse();

            while(true)
            {
               iter++;
               norm = 0;
               //F.InitF(x);
               F1.InitializationF(x);
               for(int i = 0;i < size;i++)
                F[0,i] = F1[i];
               var a = J.SolutionSystem(F1);
                for(int i = 0;i < size;i++)
                    currentX[i] = a[0,i];
                for(int i = 0;i < size;i++)
                {
                   xk1[i] = x[i] - currentX[i];

                   //norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);

                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    numOfOperations += 5;
                }
                
                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Число арифметических операций :" + numOfOperations + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");

            return xk1;
         }   
    
        // Method with transition to modified 
        public double[] TransitionToModMeth()
        {
            Console.Out.Write("Переход от обычного метода Ньютона к модифицированному :" + "\n");
            int step = 8;
            Console.WriteLine("Шаг k = " + step);
            x.InitX();

            double[] xFix = new double[size];
            xFix.InitX();

            double[] x0 = new double[size];
            x0.InitX();
            iter = 0;
            numOfOperations = 0;
            sw.Start();
            while(true)
            {
                iter++;
                norm = 0;

                F1.InitializationF(x);
                for(int i = 0;i < size;i++)
                    F[0,i] = F1[i];

                if(step > 0)
                {
                    J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);
                    J.LUDecomposition();
                    numOfOperations += J.numOfOperations;

                    for(int i = 0;i < size;i++)
                        xFix[i] = x[i];
                    //J.Reverse();
                    currentX = SquareMatrix.multiplyVec(J.reverseMatrix,F1,size);
                    step--;
                }
                else
                {
                    if(step == 0)
                    {
                        J = new SquareMatrix(InitializationExtension.InitializationJ(xFix),size);
                        J.LUDecomposition();
                    }
                    numOfOperations += J.numOfOperations;
                    //currentX = J.SolutionSystem(F);
                    currentX = SquareMatrix.multiplyVec(J.reverseMatrix,F1,size);

                    step--;
                }
                for(int i = 0;i < size;i++)
                {
                    xk1[i] = x[i] - currentX[i];

                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    numOfOperations += 5;
                }
                if(norm < eps)
                    break;

            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Число арифметических операций :" + numOfOperations + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");
            return xk1;
        }

        // method with recounting reverse matrix every k iterations
        public double[] MethodWithRewriteRevMat()
        {
            int step = 7;
            Console.Out.Write("Поиск обратной матрицы каждые k = " + step + " итераций :" + "\n");
            
            
            x.InitX();
            //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
            J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);

            J.LUDecomposition();
           numOfOperations += J.numOfOperations;

            iter = 0;
            numOfOperations = 0;
            sw.Start();
            while(true)
            {
                iter++;
                norm = 0;
                F1.InitializationF(x);
                for(int i = 0;i < size;i++)
                    F[0,i] = F1[i];
                if(iter % step == 1)
                {
                    //J = new SquareMatrix(InitializationExtension.InitJ(x),size);
                J = new SquareMatrix(InitializationExtension.InitializationJ(x),size);

                    J.LUDecomposition();
                    numOfOperations += J.numOfOperations;

                    J.Reverse();
                }
                var a = J.SolutionSystem(F1);
                for(int i =0;i < size;i++)
                    currentX[i] = a[0,i];

                for(int i = 0;i < size;i++)
                {
                    xk1[i] = x[i] - currentX[i];

                    //norm += (xk1[i] - x[i]) * (xk1[i] - x[i]);

                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    numOfOperations += 5;
                }

                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Количество итераций :" + iter + "\n");
            Console.Out.Write("Число арифметических операций :" + numOfOperations + "\n");
            Console.Out.Write("Затраченное время :" + sw.Elapsed + "\n");
            return xk1;
        }
    }
}