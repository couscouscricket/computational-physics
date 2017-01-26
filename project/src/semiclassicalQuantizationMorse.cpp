#include <iostream>
#include <cmath>
#include <sstream>

// H_2 molecule
const double GAMMA = 21.7;
const double V0 = 4.747; // eV
const double r_min = 0.74166; // Angstroms
const double E0 = -4.477;

inline double V(double r, double beta)
{
  // Morse potential.
  double f = 1 - exp((r_min-r)/beta);
  return V0 * ( f*f - 1);

}
double action(double E, double beta)
{
  // get turning points analytically for lennard-jones.
  double r_in = r_min - beta*log(+sqrt(E/V0+1) + 1);
  double r_out = r_min - beta*log(-sqrt(E/V0+1) + 1);

  int N = 128;
  double h = (r_out-r_in)/N;
  // integrate using bode's rule
  double y_h = 0; // E - V(r_in,beta) = 0
  for (int j = 1; j < N; j+=2)
    y_h += 32*sqrt(E - V(r_in + j*h, beta));
  for (int j = 2; j < N; j+=4)
    y_h += 12*sqrt(E - V(r_in + j*h, beta));
  for (int j = 4; j < N; j+=4)
    y_h += 14*sqrt(E - V(r_in + j*h, beta));
  return y_h *= GAMMA*2*h/45;
}

// action is equal to de broglie wave number,
// so define funct and find its roots.
inline double funct(double E, double beta, int n){ return action(E,beta) - (n + .5)*2*M_PI; }
double energy(double beta, int n)
{
  const double step = 0.001; //2^-9
  double E1 = E0;
  double x = funct(E1,beta,n);
  bool flag = (x>0);
  for (double E2 = E1+step; E2 < 0; E2+=step)
  {
    x = funct(E2,beta,n);
    if ((x>0) != flag)// function changes sign in [E1,E2]
    {
      // secant method.
      for (int i = 0; i < 10; ++i)
      {
        double y = funct(E2,beta,n);
        E2 += -y*(E2-E1)/(y-funct(E1,beta,n));
      }
      return E2;
    }
    E1 = E2;
  }
  return 1.;
}

int main()
{
  const double beta = .73318636;
  std::ostringstream str_gp;
  str_gp << "set terminal epslatex standalone\n";
  str_gp << "set output 'thisWillBeErased.tex'\n";
  str_gp << "set colorsequence podo\n";
  str_gp << "set border lw 3\n";
  str_gp << "set sample 300\n";
  str_gp << "set key top right\n";
  // str_gp << "set xtics 0,1,14\n";
  // str_gp << "set xlabel '$n$'\n";
  // str_gp << "set ylabel '$E_n$'\n";
  // str_gp << "set xrange [-1:15]\n";
  // str_gp << "set yrange [0:-5]\n";
  str_gp << "plot ";

  // plots action
  str_gp << "'-' w l lw 3 t '$S(E)$'\n";

  // // plot quantized energies.
  // str_gp << "'-'  w p pt 7 ps 1.5 t '$\\textbf{roots: }S(E_n)-(n+1/2)2\\pi$',";
  // str_gp << "'-'  w labels notitle\n";

  ///////////////////////////////////////////////////////////////////
  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "%s", str_gp.str().c_str());

  // plots action.
  for (double E = -V0; E < 0; E+=0.001953125)
    fprintf(gp, "%f %f\n", E, action(E,beta));
  fprintf(gp, "e\n");
  
  // plots quantized energies.
  // for (int n = 0; n < 15; ++n)
  //   fprintf(gp, "%i %f\n", n, energy(beta,n));
  // fprintf(gp, "e\n");
  // for (int n = 0; n < 15; ++n)
  // {
  //   double E = energy(beta,n);
  //   fprintf(gp, "%i %f %.3f\n", n, E-.2, E);
  // }
  // fprintf(gp, "e\n");

  pclose(gp);
  // /////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////
  std::string str_sys = "";
  str_sys += "latex -interaction batchmode thisWillBeErased.tex\n";
  str_sys += "dvipdf thisWillBeErased.dvi output.pdf\n";
  str_sys += "rm -f thisWillBeErased*\n";
  system(str_sys.c_str());
  return 0;
}