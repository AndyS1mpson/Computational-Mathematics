using System;

namespace NewtonMethod
{
    class Program
    {
        static void OutputVec(double[] x)
        {
            for (int i = 0; i < x.Length; i++)
                System.Console.Out.Write(Math.Round(x[i], 5) + " ");
            System.Console.Out.Write("\n" + "\n");
        }
        static void Main(string[] args)
        {

            NewtonMethod newton = new NewtonMethod();
            double[] result;
            result = newton.Method();
            OutputVec(result);

            result = newton.ModifiedMethod();
            OutputVec(result);

            result = newton.TransitionToModMeth();
            OutputVec(result);
            
            result = newton.MethodWithRewriteRevMat();
            OutputVec(result);
        }
    }
}
