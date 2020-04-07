namespace NewtonMethod
{
    public static class SwapExtension
    {
        public static void SwapCols(this double[,] matrix,int col1,int col2)
        {
            for(int i = 0;i < matrix.GetLength(0); i++)
            {
                double temp = matrix[i,col1];
                matrix[i,col1] = matrix[i,col2];
                matrix[i,col2] = temp;
            }
        }
        public static void SwapRows(this double[,] matrix,int str1,int str2)
        {
            for(int i = 0;i < matrix.GetLength(1);i++)
            {
                double temp = matrix[str1,i];
                matrix[str1,i] = matrix[str2,i];
                matrix[str2,i] = temp;
            }
        }
    }
}