namespace SolvingSLAE
{
    public static class ReverseExtension
    {
        public static double[,] Transpose(this double[,] matrix,int size)
        {
            double[,] transpMatrix = new double[size,size];
            for(int i = 0;i < size;i++)
            {    
                transpMatrix[i,i] = matrix[i,i];
                for(int j = 0;j < i;j++)
                {
                    transpMatrix[i,j] = matrix[j,i];
                    transpMatrix[j,i] = matrix[i,j];
                }
            }
        return transpMatrix;
        }
    }
}