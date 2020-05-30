using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Numerics;


namespace DifferentialEquation
{
    class DrawGraph
    {
        [DllImport("user32.dll")]
        public static extern IntPtr GetDC(IntPtr hWnd);

        [DllImport("kernel32.dll")]
        public static extern IntPtr GetConsoleWindow();

        [DllImport("user32.dll")]
        public static extern bool ReleaseDC(IntPtr hWnd, IntPtr hDc);

        [DllImport("gdi32.dll")]
        static extern IntPtr DeleteDC(IntPtr hDc);

        static IntPtr hWnd = IntPtr.Zero;
        static IntPtr hDC = IntPtr.Zero;

        public static void Draw(int xStart, int yStart, Complex[,] grid1, Complex[,] mas2, double[] mist)
        {
            var b = new Bitmap(1500, 2000);
            using (Graphics g = Graphics.FromImage(b))
            {

                Pen blackPen = new Pen(Color.Black);
                Pen brownPen = new Pen(Color.Brown, 2);
                Pen greenPen = new Pen(Color.DarkGreen, 2);
                Font font = new Font("verdana", 13);

                g.FillRectangle(Brushes.White, 0, 0, 1500, 2000);
                g.DrawLine(blackPen, new Point(20, 70), new Point(300, 70));
                g.DrawLine(blackPen, new Point(30, 70), new Point(30, 300));
                g.DrawString("X", font, Brushes.White, 310, 10);
                g.DrawString("Y", font, Brushes.White, 40, 250);


                int y = 3;
                int x = 100;
                for (; x < 1300; x++)
                {
                    y = (int)(grid1[0, x].Real * 70 + 70);
                    g.DrawRectangle(brownPen, x, y, 1, 1);
                }

                g.DrawLine(blackPen, new Point(20, 500), new Point(300, 500));
                g.DrawLine(blackPen, new Point(30, 500), new Point(30, 800));
                g.DrawString("X", font, Brushes.White, 310, 10);
                g.DrawString("Y", font, Brushes.White, 40, 250);
                x = 100;
                for (; x < 1300; x++)
                {
                    y = (int)(mas2[0, x].Real * 70 + 500);
                    g.DrawRectangle(greenPen, x, y, 1, 1);
                }

                g.DrawLine(blackPen, new Point(20, 1000), new Point(300, 1000));
                g.DrawLine(blackPen, new Point(30, 1000), new Point(30, 1300));
                g.DrawString("X", font, Brushes.Black, 310, 990);
                g.DrawString("Y", font, Brushes.Black, 40, 1250);
                x = 0;
                for (; x < 6; x++)
                {
                    y = (int)(mist[x] / 10 + 1000);
                    g.DrawRectangle(greenPen, x * 20 + 100, y, 5, 5);
                }


                font.Dispose();
                blackPen.Dispose();
                brownPen.Dispose();
            }
            b.Save("Graphs.png");

        }

    }
}
