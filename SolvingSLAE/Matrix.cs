using System;

namespace SolvingSLAE
{
    public class Matrix
    {
        public double eps;
        public double[,] _matrix;              // main matrix A
        public double[,] matrixL,matrixU;      // matrices of LU decomposition
        public double[,] matrixP,matrixQ;      // matrices : LU = PAQ
        public int sizeCols,sizeRows;          // size of matrix A
        public bool isLU;                       // the check the LU decomposition exist or not
        public bool isSingular;                 // the check matrix is singular or not
        private int _rank;                      // the rank of matrix
        public double norm;
        public Matrix(double[,] matrix)
        {

            sizeCols = matrix.GetLength(1);
            sizeRows = matrix.GetLength(0);

            _matrix = new double[sizeRows,sizeCols];

            for(int i = 0;i < sizeRows;i++)
            for(int j = 0;j < sizeCols;j++)
                _matrix[i,j] = matrix[i,j];
            isLU = false;

            matrixL = new double[sizeRows,sizeRows];
            matrixU = new double[sizeRows,sizeCols];
            matrixP = new double[sizeRows,sizeRows];
            matrixQ = new double[sizeCols,sizeCols];
            
            if(sizeCols > sizeRows)
                _rank = sizeRows;
            else _rank = sizeCols;

            for(int i =0;i < sizeRows;i++)
                matrixP[i,i] = 1;
            for(int j = 0;j < sizeCols;j++)
                matrixQ[j,j] = 1; 
            norm = Norm(_matrix);           
            eps = 1e-14 * norm;
        }  

        public int GetRank()
        {
            if(isLU)
                return _rank;
            else
            {
                LUDecomposition();
                return _rank;
            }
        }
           public double Norm(double[,] matrix)
        {
            double cur;
            double _norm = 0;
            for(int i = 0;i < matrix.GetLength(0);i++)
            {
                cur = 0;
                for(int j = 0;j < matrix.GetLength(1);j++)
                {
                    cur += Math.Abs(matrix[i,j]);
                }
                if(cur > norm)
                    _norm = cur;
            }
            return _norm;
        } 
        public void LUDecomposition()
        {
            isLU = true;
            isSingular = true;
            for(int i = 0;i < sizeRows;i++)
            {
                for(int j = 0;j < sizeCols;j++)
                {
                    matrixU[i,j] = _matrix[i,j];
                }
            }

            double leadEl = 0;
            int lineNum1,lineNum2;
            for(int i = 0;i < sizeRows;i++)
            {
                if(i < sizeCols)
                {
                    if(matrixU[i,i] == 0)
                    {
                        isSingular = false;
                        if(i == sizeRows - 1)
                        {
                            for(int j = i + 1;j < sizeCols;j++)
                            {
                                if(matrixU[i,j] != 0)
                                {
                                    isSingular = true;
                                    break;
                                }
                            }
                        }
                        for(int j = i + 1;j < sizeRows;j++)
                        {
                            if(matrixU[j,i] != 0)
                            {
                                lineNum1 = i;
                                lineNum2 = j;
                                matrixU.SwapRows(lineNum1,lineNum2);
                                matrixL.SwapRows(lineNum1,lineNum2);
                                matrixP.SwapRows(lineNum1,lineNum2);
                                i = i - 1;
                                isSingular = true;
                                break;
                            }
                        }
                        if(i >= _rank)
                            break;
                        if(isSingular == false)
                        {
                            if(_rank - 1 > i)
                            {
                                matrixU.SwapCols(i,_rank - 1);
                                matrixQ.SwapCols(i,_rank - 1);
                            }
                            _rank--;
                            i = i - 1;
                        }
                    }
                    else
                    {
                    leadEl = matrixU[i,i];
                    lineNum1 = i;
                    lineNum2 = i;
                    for(int j = i;j < sizeRows;j++)
                    {
                        if(Math.Abs(matrixU[j,i]) > leadEl && matrixU[j,i] != 0)
                        {
                            leadEl = matrixU[j,i];
                            lineNum2 = j;
                        }

                    }
                    if(lineNum1 != lineNum2)
                    {
                        matrixU.SwapRows(lineNum1,lineNum2);
                        matrixP.SwapRows(lineNum1,lineNum2);
                        matrixL.SwapRows(lineNum1,lineNum2);
                    }

                    double coeff;
                    for(int k = i + 1;k < sizeRows;k++)
                    {
                        if(matrixU[i,i] != 0)
                        {
                            coeff = matrixU[k,i] / matrixU[i,i];
                            for(int j = i;j < sizeCols;j++)
                            {
                                matrixU[k,j] -= matrixU[i,j] * coeff;
                                if(Math.Abs(matrixU[k,j]) < eps)
                                    matrixU[k,j] = 0;
                            }
                            matrixL[k,i] = coeff;
                        }
                    }
                    }
                }
            }
                for(int i = 0;i < sizeRows;i++)
                    matrixL[i,i] = 1;
        }

        public static void Output(double[,] matrix,int sizeRows,int sizeCols)
        {
            for(int i = 0;i < sizeRows;i++)
            {
                for(int j = 0;j < sizeCols;j++)
                    Console.Out.Write(Math.Round(matrix[i,j],5) + " ");
                Console.Out.Write("\n");
            }
            Console.Out.Write("\n" + "\n");
        }

        public static double[,] Multiply(double[,] matr1,double[,] matr2,int sizeRows1,int c,int sizeCols2)
        {
            double temp = 0;
            double[,] matrix = new double[sizeRows1,sizeCols2];
            for(int i = 0;i < sizeRows1;i++)
            {
                for(int j = 0;j < sizeCols2;j++)
                {
                    for(int k = 0;k < c;k++)
                    {
                        temp += matr1[i,k] * matr2[k,j];
                    }
                    matrix[i,j] = temp;
                    temp = 0;
                }
            }
            return matrix;
        }
        
        public double[] multiplyVec(double[,] matrix,double[] vec,int vecSize)
        {
            double[] result = new double[matrix.GetLength(0)];
            for(int i = 0;i < matrix.GetLength(0);i++)
            {
                for(int j = 0;j < matrix.GetLength(1);j++)
                {
                     result[i] += vec[j] *   matrix[i,j];
                }
            }
            return result;
        }
        
        public double[,] Solution(double[] b)
        {
            double[,] resultSolution;
            if(!isLU)
                LUDecomposition();
            int minSize;
            if(sizeCols > sizeRows)
            {
                minSize = sizeRows;
            }
            else 
            {
                minSize = sizeCols;
            }
            double[] X = new double[sizeCols];
            double[] Y = new double[sizeRows];
 
                int fsrNum = sizeCols - _rank;
                if(fsrNum !=0)
                resultSolution = new double[fsrNum,sizeCols];
                else{
                resultSolution = new double[1,sizeCols];
                fsrNum = 1;
                }
                
                Random rnd = new Random();
                for(int i = 0;i < fsrNum;i++)
                {
                    for(int j = sizeCols - 1;j >_rank - 1;j--)
                    {
                        X[j] = rnd.Next(0,2);
                    }
                    // look for y :
                    // Ly = b
                    for(int j = 0;j < sizeRows;j++)
                    {
                        Y[j] = b[j];
                        for(int k = 0;k < j;k++)
                            Y[j] -= Y[k] * matrixL[j,k];
                    }
                    for(int j = _rank - 1;j > -1;j--)
                    {
                        X[j] = Y[j];
                        for(int k = j + 1;k < sizeCols;k++)
                        {
                            X[j] -= X[k] * matrixU[j,k];
                        }
                        X[j] /= matrixU[j,j];
                    }
                    X = multiplyVec(matrixQ,X,minSize);
                    for(int j = 0;j < sizeCols;j++)
                        resultSolution[i,j] = X[j];
                     bool isSame = true;
                    if(i>0)
                    for(int k = i;k >0; k--)
                    for(int j = 0;j < sizeCols;j++)
                    {
                        if(resultSolution[i,j] != resultSolution[k-1,j])
                            isSame = false;
                    }
                    
                    if(isSame && i > 0)
                        i--;
                }

            return resultSolution;
        }    
    }
}