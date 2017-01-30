#include <cmath>
#include <array>
#include <vector>
#include <sstream>


typedef double function(double x, void *params);


// dimensionless quantum harmonic oscillator
// Y'' - (k^2)Y = 0
double qho(double x, void *params)
{
  double energy = *(double *) params;
  return 2*energy - x*x;
}

typedef struct numerov_params
{
  int N;
  double x_0, h, energy, *psi;
} numerov_params;
void numerov(function f, void *params)
{
  numerov_params *p = (numerov_params *) params;
  int N = p->N;
  double *psi = p->psi;
  double *e = &(p->energy);
  double *k = new double[N]();

  double h = p->h;

  double x_n = p->x_0;
  k[0] = f(x_n,e);
  k[1] = f(x_n + h,e);
  x_n += 2*h;
  double aux = 1./12 * h*h;
  for (int n = 2; n < N; ++n, x_n+=h)
  {
    k[n] = f(x_n,e);
    double a = 2.*(1. - 5 * aux * k[n-1])*psi[n-1];
    double b = (1. + aux * k[n-2])*psi[n-2];
    double c = 1. + aux * k[n];
    psi[n] = (a - b)/c;
  }
  delete[] k;
}

int main()
{
  const int N = 16384; // 2^14
  double x_i, x_f, h;
  double psi[N] = {0};

  numerov_params p = {.psi = psi, .N = N};
  for (double e = 0.; e < 100; e+=.125)
  {
    x_i = -sqrt(2*e)-1;
    x_f = +sqrt(2*e)+1;
    h = (x_f - x_i)/N;
    p.psi[1] = h;
    p.h = h;
    p.x_0 = x_i;
    p.energy = e;
    numerov(qho,&p);
    if (fabs(psi[N-1]) < 1)
      printf("%f\n", e);
  }




  std::ostringstream gpcmd;
  gpcmd << "set terminal epslatex standalone\n";
  gpcmd << "set output 'thisWillBeErased.tex'\n";
  gpcmd << "set colorsequence podo\n";
  gpcmd << "set border lw 3\n";
  gpcmd << "set sample 300\n";
  gpcmd << "set key top left\n";
  gpcmd << "plot ";
  gpcmd << ".5 * x**2 w l lw 3 t 'quantum harmonic oscillator',";
  gpcmd << "'-' w l lw 3 lt 3 t '$\\epsilon_n$'\n";
  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "%s",gpcmd.str().c_str());//this sends all previous commands
  for (double e = 0.5; e < 100.5; e+=20.)
  {
    p.energy = e;
    x_i = -sqrt(2*p.energy)-1;
    x_f = +sqrt(2*p.energy)+1;
    h = (x_f - x_i)/N;
    p.h = h;
    p.psi[1] = h;
    p.x_0 = x_i;
    numerov(qho,&p);
    for (int n = 0; n < N; ++n,x_i+=h)
      fprintf(gp, "%f %f\n", x_i, psi[n]+p.energy);
    fprintf(gp, "\n\n");
  }
  fprintf(gp, "e\n");
  pclose(gp);

  std::string str_sys = "";
  str_sys += "latex -interaction batchmode thisWillBeErased.tex\n";
  str_sys += "dvipdf thisWillBeErased.dvi output.pdf\n";
  str_sys += "rm -f thisWillBeErased*\n";
  system(str_sys.c_str());
  return 0;
}
