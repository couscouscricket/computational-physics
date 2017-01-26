#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

inline double integrand_LJ (double x, void *params)
{
  double energy = *(double *) params;
  return sqrt(energy - 4*(gsl_pow_int(x,-12)-gsl_pow_int(x,-6)));
}
typedef struct action_params
{
  double gamma;
  gsl_function *F;
  gsl_integration_workspace *w;
} action_params;
double action(double energy, void *params)
{
  action_params *p = (action_params *) params;
  p->F->params = &energy;
  double result, error;
  double x_in = sqrt(cbrt( 2/energy * (+sqrt(1+energy)-1) ));
  double x_out = sqrt(cbrt( 2/energy * (-sqrt(1+energy)-1) ));
  gsl_integration_qag (p->F, x_in, x_out, 
      0, 1e-7, // epsabs, epsrel
      1000,1, //max subintervals, adaptative key(1 to 6)
      p->w, &result, &error); 
  return p->gamma * result;
}

typedef struct quantize_params
{
  int n;
  action_params *aparams;
} quantize_params;
double quantize (double e, void *params)
{
  quantize_params *p = (quantize_params *) params;
  return action(e,p->aparams) - (p->n + 0.5) * M_PI;
}
typedef struct norm_params
{
  int max_iter;
  double x_lo, x_hi;
  gsl_function *F;
  gsl_root_fsolver *s;
  quantize_params *qparams;
} norm_params;
double norm_energy(int n, void *params)
{
  norm_params *p = (norm_params *) params;
  p->qparams->n = n;
  p->F->function = &quantize;
  p->F->params = p->qparams;
  gsl_root_fsolver_set(p->s,p->F,p->x_lo,p->x_hi);
  int status, iter = 0;
  double r = 0.;
  do{
    iter++;
    status = gsl_root_fsolver_iterate (p->s);
    r = gsl_root_fsolver_root (p->s);
    p->x_lo = gsl_root_fsolver_x_lower (p->s);
    p->x_hi = gsl_root_fsolver_x_upper (p->s);
    status = gsl_root_test_interval (p->x_lo, p->x_hi,
     0, 0.001);
  }while( status == GSL_CONTINUE && iter < p->max_iter);
  return r;
}

int main ()
{
  gsl_integration_workspace *w 
  = gsl_integration_workspace_alloc (1000);

  gsl_root_fsolver *s
  = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

  gsl_function F = {.function=integrand_LJ, .params=NULL};
  action_params aparams = {.gamma=21.7, .F=&F, .w=w};
  double step = gsl_pow_int(2.,-9);
  int flag = 1; double bound = -step;
  printf("ACTION s(e):\n");
  for (double energy = -1; energy < 0; energy+=step)
  {
    printf("%f %f\n", energy, action(energy,&aparams));
    if (flag && energy >= bound)
    {
      flag = !flag;
      step = gsl_pow_int(2.,-20);
    }
  }


  gsl_function G;
  G.function = quantize;
  norm_params nparams;
  quantize_params qparams;
  nparams.max_iter = 100;
  nparams.F = &G;
  nparams.s = s;
  nparams.qparams = &qparams;
  nparams.qparams->aparams = &aparams;
  printf("\nQUANTIZED ENERGIES E_n:\n");
  for (int n = 0; n < 6; ++n)
  {
    nparams.x_lo = -0.9999;
    nparams.x_hi = -0.0001;
    printf("%i %f\n", n, 4.747*norm_energy(n,&nparams));
  }


  gsl_root_fsolver_free(s);
  gsl_integration_workspace_free (w);
  return 0;
}