using System;

namespace SolvingSLAE
{
    public class SquareMatrix
    {
        private Random rnd = new Random();
        public double eps;
        public double[,] _matrix;                   // initial matrix
        public double[,] reverseMatrix;             // inverse matrix
        public int size;                            // size of initial matrix
        public double conditionNumber;              // the condition number of matrix
        public bool isLU;                           // the check the LU decomposition exist or not
        public bool isQR;                           // the check the QR decomposition exist or not
        private bool isSingular;                    // the check if matrix is singular or not
        private double _determinant;                // the determinant of the matrix
        private int _rank;                          // the rank of the matrix
        public double[,] matrixL, matrixU;          // LU decomposition matrices
        public double[,] Q, R;                      // QR decomposition matrices
        public double[,] matrixP, matrixQ;          // unit permutation matrices
        public double norm;                         // norm of matrix
        public int jacobyInter;                     // amount of iteration in the Jacoby's method
        public int seidelIter;                      // amount of iteration in the Seidel's method
        public SquareMatrix(double[,] matrix,int size)
        {
            this.size = size;
            _matrix = matrix;
            isLU = false;
            isQR = false;
            _rank = size;
            norm = Norm(_matrix);
            eps = 1e-14;
            Reverse();
            conditionNumber = Norm(reverseMatrix) * Norm(_matrix);

        }

        // display a matrix
        public static void Output(double[,] matrix,int size)
        {
            for(int i = 0;i < size;i++)
            {
                for(int j = 0;j < size;j++)
                    Console.Out.Write(Math.Round(matrix[i,j],5) + " ");
                Console.Out.Write("\n");
            }
            Console.Out.Write("\n" + "\n");
        }

        // get determinant
        public double GetDeterminant()
        {
            if(isLU)
                return _determinant;
            else 
            {
                LUDecomposition();
                return _determinant;
            }
        }
        
        // get rank
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
        
        // look for norm of matrix
        public double Norm(double[,] matrix)
        {
            double cur;
            double _norm = 0;
            for(int i = 0;i < size;i++)
            {
                cur = 0;
                for(int j = 0;j < size;j++)
                {
                    cur += Math.Abs(matrix[i,j]);
                }
                if(cur > _norm)
                    _norm = cur;
            }
            return _norm;
        } 
        
        // implement LU matrix decomposition
        public void LUDecomposition()
        {
            isLU = true;
            isSingular = true;
            matrixL = new double[size, size];
            matrixU = new double[size, size];
            matrixP = new double[size, size];
            matrixQ = new double[size, size];

            for(int i = 0;i < size;i++)
            {
                // fill U matrix with the elements of the initial matrix
                for(int j = 0;j < size;j++)
                {
                    matrixU[i,j] = _matrix[i,j];
                }
                // make P,Q matrices singular
                matrixP[i,i] = 1;
                matrixQ[i,i] = 1;
            }

            int lineNum1,lineNum2;
            double leadElem;
            int numSwap = 0;

            for(int i = 0;i < size;i++)
            {
                // if diagonal element is 0 and next line element is not 0  make swap lines 
                if(matrixU[i,i] == 0)
                {   
                    isSingular = false;
                    for(int j = i + 1;j < size;j++)
                    {
                        if(matrixU[j,i] != 0)
                        {
                            lineNum1 = i;
                            lineNum2 = j;
                            matrixU.SwapRows(lineNum1,lineNum2);
                            matrixL.SwapRows(lineNum1,lineNum2);
                            matrixP.SwapRows(lineNum1,lineNum2);
                            i = i - 1;
                            numSwap++;
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
                            numSwap++;
                        }
                        _rank--;
                        i = i - 1;
                    }
                }
                else
                {
                    isSingular = true;
                    leadElem = matrixU[i,i];
                    lineNum1 = i;
                    lineNum2 = i;
                    // look for the max element in a column
                    for(int j = i;j < size;j++)
                    {
                        if(Math.Abs(matrixU[j,i]) > leadElem && matrixU[j,i] != 0)
                        {
                            leadElem = Math.Abs(matrixU[j,i]);
                            lineNum2 = j;
                        }
                    }
                    // swap the line with leadElem and current line
                    if(lineNum1 != lineNum2)
                    {
                        matrixU.SwapRows(lineNum1,lineNum2);
                        matrixL.SwapRows(lineNum1,lineNum2);
                        matrixP.SwapRows(lineNum1,lineNum2);
                        numSwap++;
                    }
                    double coeff;        // the coeff for zeroing column's elements under current diagonal element
                    for(int k = i + 1;k < size;k++)
                    {
                        if(matrixU[i,i] != 0)
                        {
                            coeff = matrixU[k, i] / matrixU[i, i];
                            for(int j = i;j < size;j++)
                            {
                                matrixU[k,j] = matrixU[k, j] - matrixU[i, j] * coeff;
                                // check element is abs 0 or not
                                if(Math.Abs(matrixU[k,j]) < eps)
                                    matrixU[k,j] = 0;
                            }
                            matrixL[k, i] = coeff;
                        }
                    }
                    
                }
            }
            for(int i = 0;i < size;i++)
                matrixL[i,i] = 1;
                _determinant = 1;
            for(int i = 0;i < size;i++)
                _determinant *= matrixU[i,i];
            if(numSwap % 2 != 0)
                _determinant *= (-1);
        }
        
        // multiply matrix
        public static double[,] multiply(double[,] matrix1,double[,] matrix2,int size)
        {
            double[,] mult = new double[size,size];
            for(int i = 0;i < size;i++)
            for(int j = 0;j < size;j++)
            {
                for(int k = 0;k < size;k++)
                {
                    mult[i,j] += matrix1[i,k] * matrix2[k,j];
                }
            }
            return mult;
        }
       
        // check LU is equal to PAQ where A is _matrix
        public (double[,],double[,]) InfelicityCorrection(double[,] matrix1,double[,] matrix2)
        {
            //double[,] multLU = multiply(matrixL,matrixU,size);
            //double[,] multPAQ = multiply(multiply(matrixP,_matrix,size),matrixQ,size);
            double[,] multLU = matrix1;
            double[,] multPAQ =matrix2;
            for(int i = 0;i < size;i++)
            for(int j = 0;j < size;j++)
            {
                if(Math.Abs(multLU[i,j] - multPAQ[i,j]) < 1e-14)
                    multLU[i,j] = multPAQ[i,j];
            }
            return (multLU,multPAQ);
        } 
        
        // find solution of system
        public double[,] SolutionSystem(double[] vec)
        {
            double[,] resultSolution;
            if(!isLU)
                LUDecomposition();

            double[] X = new double[size];
            double[] Y = new double[size];
            int fsrNum = size - _rank;
                if(fsrNum !=0)
                resultSolution = new double[fsrNum,size];
                else{
                resultSolution = new double[1,size];
                fsrNum = 1;
                }


            vec = multiplyVec(matrixP,vec,size);        // conversion to view after rearrangement
            isSingular = true;

            for(int i = 0;i < size;i++)                 // system compatibility check
            {
                if(vec[i] != 0)
                {
                    isSingular = true;
                    break;
                }
            }
            for(int i = 0;i < fsrNum;i++)
            {
                for(int j = size - 1;j >_rank - 1;j--)
                        {
                            X[j] = rnd.Next(0,2) ;

                        }
                // look for y :
                // Ly = b
                for(int j = 0; j < size; j++)                //reverse stroke along the left matrix
                {
                    Y[j] = vec[j];
                        for (int k = 0; k < j; k++)
                            Y[j] -= matrixL[j,k] * Y[k];
                }
                // look for x :
                // Ux = y
                for(int j = _rank - 1;j > -1;j--)            //reverse stroke along the right matrix
                {
                    X[j] = Y[j] / matrixU[j,j];
                    for(int k = j+1;k < size;k++)
                        X[j] -= matrixU[j,k] * X[k] / matrixU[j,j];
                }
                X = multiplyVec(matrixQ,X,size);            //remove the effect of column permutation
                for(int j = 0;j < size;j++)
                    resultSolution[i,j] = X[j];
                bool isSame = true;
                if(i>0)
                for(int k = i;k >0; k--)
                for(int j = 0;j < size;j++)
                {
                    if(resultSolution[i,j] != resultSolution[k-1,j])
                        isSame = false;
                }
                
                if(isSame && i > 0)
                    i--;

            }
                // find matrix condition number
                double vecNorm = 0;
                double solutionNorm = 0;
                for(int i = 0;i < vec.Length - 1; i++)
                {
                    vecNorm += vec[i] * vec[i];
                    solutionNorm += X[i] * X[i];
                }
                vecNorm = Math.Sqrt(vecNorm);
                solutionNorm = Math.Sqrt(solutionNorm);
            return resultSolution;
        }
        
        // multiply matrix with vector
        public static double[] multiplyVec(double[,] matrix,double[] vec,int size)
        {
            double[] vector = new double[size];
            for(int i = 0;i < size;i++)
            {
                for(int j = 0;j < size;j++)
                {
                    vector[i] += vec[j] *   matrix[i,j];
                }
            }
            return vector;
        }
        
        // find inverse matrix
        public void Reverse()
        {
            if(!isLU)
                LUDecomposition();
            
            double[] vec = new double[size];
            reverseMatrix = new double[size,size];
            for(int k = 0;k < size;k++)
            {
                vec[k] = 1;
                vec = multiplyVec(matrixP,vec,size);

                for(int j = 0;j < size;j++)
                {
                    for(int i = j + 1;i < size;i++)
                        vec[i] -= matrixL[i,j] * vec[j];
                }
                for(int j = size - 1;j > -1;j--)
                {
                    vec[j] /= matrixU[j,j];
                    for(int i = 0;i < j;i++)
                        vec[i] -= matrixU[i,j] * vec[j];
                }
                for(int i = 0;i < size;i++)
                    reverseMatrix[i,k] = vec[i];
                for(int i = 0;i < size;i++)
                    vec[i] = 0;
            }

        }
        
        // QR decomposition
        // in this program the solution is implemented through the Givens rotation method 
        public void QRDecomposition()
        {
            isQR = true;
            Q = new double[size,size];
            R = new double[size,size];
            double[,] rotateMatr = new double[size,size];
            double sin,cos;                         // sine and cosine of the rotation angle of the matrix
            for(int i = 0;i < size;i++)
            {
                rotateMatr[i,i] = 1;
                Q[i,i] = 1;
                for(int j = 0;j < size;j++)
                    R[i,j] = _matrix[i,j];
            }
            // go through the columns
            for(int j = 0;j < size - 1;j++)
            {
                for(int i = j + 1;i < size;i++)
                {
                    if(R[i,j] != 0)
                    {
                        cos = R[j,j] / (Math.Sqrt( Math.Pow (R[j,j], 2) 
                                                                + Math.Pow (R[i,j], 2)));
                        sin = -R[i,j] / (Math.Sqrt( Math.Pow (R[j,j], 2) 
                                                                + Math.Pow (R[i,j], 2)));
                        
                        rotateMatr[j,j] = cos;    rotateMatr[j,i] = -sin;
                        rotateMatr[i,j] = sin;    rotateMatr[i,i] = cos;

                        R = multiply(rotateMatr,R,size);
                        Q = multiply(rotateMatr,Q,size);
                        
                        rotateMatr[j,j] = 1;    rotateMatr[j,i] = 0;
                        rotateMatr[i,j] = 0;    rotateMatr[i,i] = 1;

                    }
                }
            }
            for(int i = 0;i < size;i++)
            for(int j = 0;j < size;j++)
            {
                if(Math.Abs(Q[i,j]) < eps)
                    Q[i,j] = 0;
                if(Math.Abs(R[i,j]) < eps)
                    R[i,j] = 0;
            }
        }
        public double[] QRSolutionSystem(double[] b)
        {
            if(!isQR)
            {
                QRDecomposition();
            }
            double[] X = new double[size];
            b = multiplyVec(Q,b,size);
            for(int j = size - 1;j > -1;j--)
            {
                X[j] = b[j] / R[j,j];

            for(int i = j+1;i < size;i++)
                //b[i] -= R[i,j] * X[j];
                X[j] -= R[j,i] * X[i] / R[j,j];
            }

            return X;
        }
        
        // the method of Jacoby for solving the SLAE
        // in this method matrix A is sum of 3 matrices: L,D,R
        public double[] methodJacoby(double[] b,double[] x0)
        {
            double[,] matrixB = new double[size,size];          // B = D^(-1) * ( L + R)
            double[] Xk = new double[size];                     // kth system solution
            double[] Xk1 = new double[size];                    // (k+1)th system solution
            jacobyInter = 0;                                    // number of iteration

            Xk = x0;
            // make B = A
            for(int i = 0;i < size;i++)
            {
                for(int j = 0;j < size;j++)
                matrixB[i,j] = _matrix[i,j];
            }
            // make (L + R) matrix
            for(int i = 0;i < size;i++)
            {
                matrixB[i,i] = 0;                               
                for(int j = 0;j < size;j++)
                    matrixB[i,j] /= _matrix[i,i];
            }
            // make b[i] / A[i,i]
            for(int i = 0;i < size;i++)
                {
                    b[i] /= _matrix[i,i];
                    Xk[i] = b[i];
                }
            double dif;
            if(Norm(matrixB) < 0.5)
                dif = (1 - Norm(matrixB)) /Norm(matrixB) ;
            else dif = eps;
            
            do{
                double normDif = 0;
                jacobyInter++;
                for(int i = 0;i < size;i++)
                    Xk1[i] = Xk[i];
                // make Xk = D^(-1) * ( L + R) * Xk
                Xk = multiplyVec(matrixB,Xk,size);

                for(int i = 0;i < size;i++)
                {
                    Xk[i] = -Xk[i] + b[i];
                }

                for(int j = 0;j < size;j++)
                {
                    normDif += (Xk[j] - Xk1[j]) * (Xk[j] - Xk1[j]);
                }
                normDif = Math.Sqrt(normDif);

                if(normDif < dif)
                    break;

            }while(true);


            return Xk;
        }   
       
        // the method of Seidel for solving the SLAE
        public double[] methodSeidel(double[] b,double[] x0)
        {
            double[,] matrixB = new double[size,size];          // B = D^(-1) * ( L + R)
            double[] Xk = new double[size];                     // kth system solution
            double[] Xk1 = new double[size];                    // (k+1)th system solution
            seidelIter = 0;                                    // number of iteration
            double[] Xp = new double[size];

            for(int i = 0;i < size;i++)
                Xk[i] = x0[i];

            // make B = A
            for(int i = 0;i < size;i++)
            {
                for(int j = 0;j < size;j++)
                matrixB[i,j] = _matrix[i,j];
            }
            // make (L + R) matrix
            for(int i = 0;i < size;i++)
            {
                matrixB[i,i] = 0;                               
                for(int j = 0;j < size;j++)
                    matrixB[i,j] /= _matrix[i,i];
            }
            // make b[i] / A[i,i]
            for(int i = 0;i < size;i++)
                {
                    b[i] /= _matrix[i,i];
                    Xk[i] = b[i];
                }

            double dif;
            if(Norm(matrixB) < 0.5)
                dif = (1 - Norm(matrixB)) /Norm(matrixB) ;
            else dif = eps;

            do
            {
                seidelIter++;

                for(int i = 0;i < size;i++)
                    Xp[i] = Xk[i];

                double normDif = 0;
                
                for(int i = 0;i < size;i++)
                {
                    Xk1[i] = b[i];
                    for(int j = 0;j < size;j++)
                    {
                        Xk1[i] -= Xk[j] * matrixB[i,j];
                    }
                    Xk[i] = Xk1[i];
                }

                for(int j = 0;j < size;j++)
                    {
                        normDif += (Xk[j] - Xp[j]) * (Xk[j] - Xp[j]);
                    }
                    normDif = Math.Sqrt(normDif);

                if(normDif < dif)
                    break;


            }while(true);
            
                return Xk;
        }
    }
}