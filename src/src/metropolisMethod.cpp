#include <random>
#include <iostream>
#include <cmath>
#include <array>
#include <sstream>

inline double w(double x)
{
  //Normalized in the range [0,pi]

  return 1./M_PI * (sin(2*x)*sin(2*x) + cos(x)*cos(x));
  // return exp(-x*x)/sqrt(M_PI); //gaussian
}

int main()
{
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<> uniform(0.0, 1.);

  const int M = 100; //Partition of region of intergration
  const int walkers = M+1; //Number of walkers
  const double x2 = M_PI;
  const double x1 = 0;
  const double deltaM = (x2-x1)/M;
  const int steps = 30000;

  std::vector<double> X;
  for (int walker = 0; walker <= walkers; ++walker)
  {
    X.push_back(x1 + deltaM*walker);
    for (int s = 0; s < steps; ++s)
    {
      double Xn = X.back();
      double Xt = Xn + deltaM*(2*uniform(rng) - 1);
      if (w(Xt)/w(Xn) > uniform(rng) &&  x1 <= Xt && Xt <= x2)
        X.push_back(Xt);
    }
  }

  std::array <double,M> histogram = {};
  for (int i = 0; i < X.size(); ++i)
  {
    for (int j = 0; j < histogram.size(); ++j)
    {
      if (x1 + deltaM*j < X[i] && X[i] <= x1 + deltaM*(j+1))
        ++histogram[j];
    }
  }

    //CDF: discrete cumulative distribution function.
  std::array <double,M+1> CDF = {};
  const int totalPoints = X.size();
  for (int i = 1; i < CDF.size(); ++i)
    CDF[i] = CDF[i-1] + histogram[i-1]/totalPoints;

  ///////////// gnuplot's commands ////////////////////////////////
  std::ostringstream str_gp;
  str_gp << "set terminal epslatex standalone\n";
  str_gp << "set output 'thisWillBeErased.tex'\n";
  str_gp << "set colorsequence podo\n";
  str_gp << "set label '"<<M<<" walkers' at "<<.7*(x2-x1)+x1<<",.8\n";
  str_gp << "set label '"<<steps<<" steps' at "<<.7*(x2-x1)+x1<<",.74\n";
  str_gp << "set border lw 3\n";
  str_gp << "set key top left spacing 1.3\n";
  str_gp << "set xrange ["<<x1<<":"<<x2<<"]\n";
  str_gp << "set sample 300\n";
  str_gp << "plot ";
  //////////////////////////////////////////////////////////////////

  /////////////// plot ////////////////////////////////////////////////
  FILE *gp = popen("gnuplot","w");

  str_gp << "(sin(2*x)**2 + cos(x)**2)/pi w l lw 3 t '$w(x)$',"; //sum of cos and sin
  str_gp << "(-sin(4*x) +2*sin(2*x) +8*x)/(8*pi) w l lw 3 t '$\\int w(x)$',"; //sum of cos and sin CDF
  // str_gp << "exp(-x**2)/sqrt(pi) w l lw 3 t '$w(x)$',"; //gaussian
  // str_gp << ".5*(erf(x)+1) lw 2 t '$\\int w(x)$',"; //gaussian CDF
  str_gp << "'-' w boxes lw 3 t 'histogram',";
  str_gp << "'-' w steps lw 3 t 'cdf'\n";

  fprintf(gp, "%s", str_gp.str().c_str());//this sends all commands

  //this sends the points for command '-'
  for (int k = 0; k < histogram.size(); ++k)
    fprintf(gp, "%f %f\n", deltaM*k+x1, histogram[k]/(totalPoints*deltaM));
  fprintf(gp, "%f %f\n", x2, histogram[M-1]/(totalPoints*deltaM));
  fprintf(gp, "e\n");

  for (int k = 0; k < CDF.size(); ++k)
    fprintf(gp, "%f %f\n", deltaM*k+x1, CDF[k]);
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