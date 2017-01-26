#include <random>
#include <iostream>
#include <cmath>
#include <array>
#include <sstream>

inline double w(double x)
{
  //Normalized in the range [0,pi]

  // return 1./M_PI * (sin(2*x)*sin(2*x) + cos(x)*cos(x));
  return exp(-x*x)/sqrt(M_PI); //gaussian
}

int main()
{
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<> uniform(0.0, 1.);

  const int M = 50; //Partition of region of intergration
  const int walkers = M+1; //Number of walkers
  const double x2 = 3;
  const double x1 = -3;
  const double deltaM = (x2-x1)/M;

  const int N = 3; // Number of different walks
  std::array<int,N> stepsArray = {1, 500,1000};
  std::array<std::vector<double>,N> X;
  for (int i = 0; i < stepsArray.size(); ++i)
  {
    for (int walker = 0; walker < walkers; ++walker)
    {
      X[i].push_back(x1 + deltaM*walker);
      for (int s = 1; s < stepsArray[i]; ++s)
      {
        double Xn = X[i].back();
        double Xt = Xn + deltaM*(2*uniform(rng) - 1);
        if (w(Xt)/w(Xn) > uniform(rng) &&  x1 <= Xt && Xt <= x2)
          X[i].push_back(Xt);
      }
    }
  }

  std::array<std::array<double,M>,N> histogram = {};
  for (int i = 0; i < X.size(); ++i)
  {
    for (int j = 0; j < X[i].size(); ++j)
    {
      for (int k = 0; k < histogram[i].size(); ++k)
      {
        if (x1 + deltaM*k < X[i][j] && X[i][j] <= x1 + deltaM*(k+1))
          ++histogram[i][k];
      }
    }
  }

    //CDF: discrete cumulative distribution function.
  std::array<std::array<double,M+1>,N> CDF = {};
  for (int i = 0; i < X.size(); ++i)
  {
    for (int j = 1; j < CDF[i].size(); ++j)
      CDF[i][j] = CDF[i][j-1] + histogram[i][j-1]/X[i].size();
  }

  ///////////// gnuplot's commands ////////////////////////////////
  std::ostringstream str_gp;
  str_gp << "set terminal epslatex standalone\n";
  str_gp << "set output 'thisWillBeErased.tex'\n";
  str_gp << "set grid x y z back\n";
  str_gp << "set ytics (1000,15000,30000,60000)\n";
  str_gp << "set xyplane .1\n";
  str_gp << "set ylabel 'walk size' rotate parallel\n";
  str_gp << "set colorsequence podo\n";
  str_gp << "set border lw 3\n";
  str_gp << "set nokey\n";
  str_gp << "set xrange ["<<x1<<":"<<x2<<"]\n";
  str_gp << "set sample 300\n";
  str_gp << "splot ";
  //////////////////////////////////////////////////////////////////

  /////////////// plot ////////////////////////////////////////////////
  FILE *gp = popen("gnuplot","w");

  str_gp << "'-' w l lw 3,";
  str_gp << "'-' w l lt 3 dt 3 lw 5\n";

  fprintf(gp, "%s", str_gp.str().c_str());//this sends all commands
  for (int i = 0; i < N; ++i)
  {
    for (int k = 0; k < histogram[i].size(); ++k)
      fprintf(gp, "%f %i %f\n", deltaM*k+x1, stepsArray[i], histogram[i][k]/(X[i].size()*deltaM));
    fprintf(gp, "%f %i %f\n\n\n", x2, stepsArray[i], histogram[i][M-1]/(X[i].size()*deltaM));
  }
  fprintf(gp, "e\n");
  for (int i = 1; i < N; ++i)
  {
    for (double x=x1; x <= x2; x+=0.015625)
    {
      // double sine = sin(2*x);
      // double cosine = cos(x);
      // fprintf(gp, "%f %i %f\n", x, stepsArray[i], (sine*sine + cosine*cosine)/M_PI);
      fprintf(gp, "%f %i %f\n", x, stepsArray[i], exp(-x*x)/sqrt(M_PI));
    }
    fprintf(gp, "\n\n");
  }
  fprintf(gp, "e\n");
  pclose(gp);
  //////////////////////////////////////////////////////////////////

  //////////// tex2pdf and output cleanup ////////////////////////////
  std::string str_sys = "";
  str_sys += "latex -interaction batchmode thisWillBeErased.tex\n";
  str_sys += "dvipdf thisWillBeErased.dvi output.pdf\n";
  str_sys += "rm -f thisWillBeErased*\n";
  const char * SYSCMD = str_sys.c_str();
  system(SYSCMD);
  ///////////////////////////////////////////////////////////////////
  return 0;
}