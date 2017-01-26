#include <iostream>
#include <cmath>

// H_2 molecule
const double GAMMA = 21.7;
const double V0 = 4.747; // eV

inline double v(double x)
{
  // Lennard-Jones potential.
  return 4*(pow(x,-12)-pow(x,-6));
}
double action(double e)
{
  // get turning points analytically for lennard-jones.
  double x_in = sqrt(cbrt( 2/e * (+sqrt(1+e)-1) ));
  double x_out = sqrt(cbrt( 2/e * (-sqrt(1+e)-1) ));

  int N = 8192;
  double h = (x_out-x_in)/N;
  // integrate using bode's rule
  double y_h = 0; // e - v(x_in) = 0
  for (int j = 1; j < N; j+=2)
    y_h += 32*sqrt(e - v(x_in + j*h));
  for (int j = 2; j < N; j+=4)
    y_h += 12*sqrt(e - v(x_in + j*h));
  for (int j = 4; j < N; j+=4)
    y_h += 14*sqrt(e - v(x_in + j*h));
  return y_h *= GAMMA*2*h/45;
}

// action is equal to de broglie wave number,
// so define funct and find its roots.
inline double funct(double e, int n){ return action(e) - (n + .5)*M_PI; }
double normalized_energy(int n)
{
  const double step = 0.001953125; //2^-9
  double e1 = -1;
  double x = funct(e1,n);
  bool flag = (x>0);
  for (double e2 = e1+step; e2 < 0; e2+=step)
  {
    x = funct(e2,n);
    if ((x>0) != flag)// function changes sign in [e1,e2]
    {
      // secant method.
      for (int i = 0; i < 10; ++i)
      {
        double y = funct(e2,n);
        e2 += -y*(e2-e1)/(y-funct(e1,n));
      }
      return e2;
    }
   e1 = e2;
  }
  return 1.;
}

int main()
{
  for (double e = -1; e < 0; e+=0.001953125)
   printf ("%f %f\n", e, action(e));


  std::string str_gp = "";
  str_gp += "set terminal epslatex standalone\n";
  str_gp += "set output 'thisWillBeErased.tex'\n";
  str_gp += "set colorsequence podo\n";
  str_gp += "set border lw 3\n";
  str_gp += "set sample 300\n";
  str_gp += "set key top right\n";
  str_gp += "set xtics (0,1,2,3,4)\n";
  str_gp += "set xlabel '$n$'\n";
  str_gp += "set ylabel '$E_n = \\epsilon_n V_0$'\n";
  str_gp += "set xrange [-.5:4.5]\n";
  str_gp += "set yrange [0:-5]\n";
  str_gp += "plot ";

  // plots action
  // str_gp += "'-' w l lw 3 t '$s(\\epsilon)$'\n";

  // plot quantized energies.
  str_gp += "'-'  w p pt 7 ps 1.5 t '$\\textbf{roots: }s(\\epsilon_n)-(n+1/2)\\pi$',";
  str_gp += "'-'  w labels notitle\n";

  ///////////////////////////////////////////////////////////////////
  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "%s", str_gp.c_str());

  // plots action.
  // for (double e = -1; e < 0; e+=0.001953125)
  //   fprintf(gp, "%f %f\n", e, action(e));
  // fprintf(gp, "e\n");
  
  // plots quantized energies.
  for (int n = 0; n < 5; ++n)
    fprintf(gp, "%i %f\n", n, V0*normalized_energy(n));
  fprintf(gp, "e\n");
  for (int n = 0; n < 5; ++n)
  {
    double energy = V0*normalized_energy(n);
    fprintf(gp, "%i %f %f\n", n, energy-.2, energy);
  }
  fprintf(gp, "e\n");

  pclose(gp);
  /////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////
  std::string str_sys = "";
  str_sys += "latex -interaction batchmode thisWillBeErased.tex\n";
  str_sys += "dvipdf thisWillBeErased.dvi output.pdf\n";
  str_sys += "rm -f thisWillBeErased*\n";
  system(str_sys.c_str());
  return 0;
}