using System;

namespace SolvingSLAE
{
    class Program
    {
        public static double[,] inputMatrix(int sizeRows,int sizeCols)
        {
            double[,] matrix = new double[sizeRows,sizeCols];
            string[] str;
            for (int j = 0; j < sizeRows; j++)
            {
                str = Console.ReadLine().Split(new char[] { ' ', '\n', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                for (int i = 0; i < sizeCols; i++)
                {
                    matrix[j, i] = Double.Parse(str[i]);
                }
            }
            return matrix;
        } 
        static void Main(string[] args)
        {
            #region For Square Matrix
            //Console.WriteLine("Введите размер n матрицы A [n x n] :");
            Random rand = new Random();
            int size = rand.Next(2, 7);
            Console.WriteLine("Размерность матрицы : {0}", size);
            double[,] arr = new double[size, size];
            //Console.WriteLine("Введите матрицу A размерности n :");
            //arr = inputMatrix(size,size);

            // generation matrix with diagonal dominance
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    if (i == j)
                        arr[i, j] = rand.Next(50, 100);
                    else
                        arr[i, j] = rand.Next(-49, 49);
                }


            Console.WriteLine("Матрица A с диагональным преобладанием:");
            Matrix.Output(arr, size, size);

            SquareMatrix matrix = new SquareMatrix(arr, size);

            #region LU decomposition
            matrix.LUDecomposition();
            Console.WriteLine("Матрица L :");
            SquareMatrix.Output(matrix.matrixL, size);
            Console.WriteLine("Матрица U :");
            SquareMatrix.Output(matrix.matrixU, size);
            Console.WriteLine("Произведение матриц L и U :");
            SquareMatrix.Output(SquareMatrix.multiply(matrix.matrixL, matrix.matrixU, size), size);
            Console.WriteLine("Произведение матриц P,A и Q :");
            SquareMatrix.Output(SquareMatrix.multiply(SquareMatrix.multiply(matrix.matrixP, matrix._matrix, size), matrix.matrixQ, size), size);
            #endregion

            #region determinant, rank and reverse matrix
            // determinant
            Console.WriteLine("Определитель матрицы A :");
            Console.WriteLine(matrix.GetDeterminant());

            // rank
            Console.WriteLine("Ранг матрицы :");
            Console.WriteLine(matrix.GetRank());

            // the reverse matrix
            Console.WriteLine("Обратная матрица :");
            double[,] reverseM = matrix.reverseMatrix;
            SquareMatrix.Output(reverseM, size);

            Console.WriteLine(" A * A^(-1) :");
            double[,] mult = SquareMatrix.multiply(arr, matrix.reverseMatrix, size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (mult[i, j] < matrix.eps)
                        mult[i, j] = 0;
            SquareMatrix.Output(mult, size);


            Console.WriteLine(" A^(-1) * A :");
            mult = SquareMatrix.multiply(matrix.reverseMatrix, arr, size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (mult[i, j] < matrix.eps)
                        mult[i, j] = 0;
            SquareMatrix.Output(mult, size);

            #endregion

            #region solving the SLAE through LU decomposition and the condition number

            Console.WriteLine("Введите вектор b : Ax = b :");
            double[] b = new double[size];
            var m = inputMatrix(1, size);
            for (int i = 0; i < size; i++)
                b[i] = m[0, i];
            Console.WriteLine("Решение СЛАУ :");
            double[] x = matrix.SolutionSystem(b);
            for (int i = 0; i < size; i++)
                Console.WriteLine(x[i]);

            Console.WriteLine("Проверка Ax - b = 0 :");
            Console.WriteLine(" Ax : ");
            double[] check = SquareMatrix.multiplyVec(arr, x, size);
            for (int i = 0; i < size; i++)
                Console.WriteLine(check[i]);
            bool result = true;
            for (int i = 0; i < size; i++)
            {
                check[i] -= b[i];
                if (check[i] < matrix.eps)
                    check[i] = 0;
            }
            for (int i = 0; i < size; i++)
                if (check[i] != 0)
                    result = false;
            if (result)
                Console.WriteLine("Все верно!");
            else Console.WriteLine("Не верно!");

            // the condition number
            Console.WriteLine("Число обусловленности :");
            Console.WriteLine(matrix.conditionNumber);

            #endregion

            #region method QR decomposition

            matrix.QRDecomposition();
            Console.WriteLine("Матрица Q :");
            SquareMatrix.Output(matrix.Q, size);
            Console.WriteLine("Матрица R :");
            SquareMatrix.Output(matrix.R, size);
            Console.WriteLine("Произведение матриц Q и R :");
            SquareMatrix.Output(SquareMatrix.multiply(matrix.Q.Transpose(size), matrix.R, size), size);

            #endregion

            #region solving the SLAE through QR decomposition

            Console.WriteLine("Введите вектор b : Ax = b :");
            double[] b1 = new double[size];
            var m2 = inputMatrix(1, size);
            for (int i = 0; i < size; i++)
                b1[i] = m2[0, i];
            Console.WriteLine("Решение СЛАУ через QR разложение :");
            double[] x1 = matrix.QRSolutionSystem(b);
            for (int i = 0; i < size; i++)
                Console.WriteLine(x1[i]);

            #endregion

            #region methods Jacoby and Seildel 

            Console.WriteLine("Введите b : Ax = b :");
            double[] b2 = new double[size];
            var m3 = inputMatrix(1, size);
            for (int i = 0; i < size; i++)
                b2[i] = m3[0, i];
            Console.WriteLine("Введите начальное значение (x0) :");
            double[] x0 = new double[size];
            var mX = inputMatrix(1, size);
            for (int i = 0; i < size; i++)
                x0[i] = mX[0, i];
            double[] copB = new double[size];
            for (int i = 0; i < size; i++)
                copB[i] = b[i];
            double[] solutionJ = matrix.methodJacoby(b, x0);
            Console.WriteLine("Решение методом Якоби :");
            for (int i = 0; i < size; i++)
                Console.WriteLine(solutionJ[i]);
            Console.WriteLine("Число итераций :");
            Console.WriteLine(matrix.jacobyInter);

            double[] solutionS = matrix.methodSeidel(copB, x0);
            Console.WriteLine("Решение методом Зейделя :");
            for (int i = 0; i < size; i++)
                Console.WriteLine(solutionS[i]);
            Console.WriteLine("Число итераций :");
            Console.WriteLine(matrix.seidelIter);

            #endregion
            #endregion

            #region For a matrix of any dimension

            ////Console.WriteLine("Введите количество строк :");
            ////int sizeRows = Convert.ToInt32(Console.ReadLine());
            ////Console.WriteLine("Введите количество столбцов :");
            ////int sizeCols = Convert.ToInt32(Console.ReadLine());
            //Random rand = new Random();
            //int sizeRows = rand.Next(2,10);
            //int sizeCols = rand.Next(2,10);

            //double[,] arr2 = new double[sizeRows, sizeCols];
            ////Console.WriteLine("Введите матрицу A : ");
            ////arr2 = inputMatrix(sizeRows, sizeCols);

            //// matrix generation 
            //for (int i = 0; i < sizeRows; i++)
            //    for (int j = 0; j < sizeCols; j++)
            //        arr2[i, j] = rand.Next(-100, 100);




            //Matrix.Output(arr2,sizeRows,sizeCols);
            //Matrix matr1 = new Matrix(arr2);



            //// LU decomposition
            //matr1.LUDecomposition();
            //Console.WriteLine("Матрица L :");
            //Matrix.Output(matr1.matrixL, matr1.matrixL.GetLength(0), matr1.matrixL.GetLength(1));
            //Console.WriteLine("Матрица U :");
            //Matrix.Output(matr1.matrixU, matr1.matrixU.GetLength(0), matr1.matrixU.GetLength(1));
            //Console.WriteLine("Произведение матриц L и U :");
            //Matrix.Output(Matrix.Multiply(matr1.matrixL, matr1.matrixU, matr1.matrixL.GetLength(0)
            //                               , matr1.matrixL.GetLength(1), matr1.matrixU.GetLength(1)),
            //                               matr1.matrixL.GetLength(0), matr1.matrixU.GetLength(1));
            //Console.WriteLine("Произведение матриц P,A и Q :");
            //Matrix.Output(Matrix.Multiply(Matrix.Multiply(matr1.matrixP, matr1._matrix, matr1.matrixP.GetLength(0)
            //                                           , matr1.matrixP.GetLength(1), matr1._matrix.GetLength(1))
            //                                           , matr1.matrixQ, matr1.sizeRows, matr1.matrixQ.GetLength(0), matr1.matrixQ.GetLength(1)), matr1.sizeRows, matr1.sizeCols);
            //// rank
            //Console.WriteLine("Ранг матрицы A :");
            //Console.WriteLine(matr1.GetRank());

            //// solution of the SLAE
            //Console.WriteLine("Введите b : Ax = b :");
            //double[] b3 = new double[sizeRows];
            //var m4 = inputMatrix(1, sizeRows);
            //for (int i = 0; i < sizeRows; i++)
            //    b3[i] = m4[0, i];
            //double[,] sol = matr1.Solution(b3);
            //Console.WriteLine("Решение системы :");
            //for (int i = 0; i < sol.GetLength(0); i++)
            //    for (int j = 0; j < sol.GetLength(1); j++)
            //        Console.WriteLine(sol[i, j]);

            #endregion

        }
    }
}
