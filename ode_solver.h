//solver for ode Eqs.
//copied from https://people.sc.fsu.edu/~jburkardt/c_src/rk4/rk4.html
//Thanks author John Burkardt
#ifndef __ODE_SOLVER_H
#define __ODE_SOLVER_H
#include <stdlib.h>

int rk4vec ( double t0, int n, const double * u0, double * u, double dt, 
  double *f ( double t, const double * u ) );
#endif /* #ifndef __ODE_SOLVER_H */
