#include "ode_solver.h"

/******************************************************************************/

int rk4vec ( double t0, int m, const double * u0, double * u, double dt, 
  double *f ( double t, const double * u ) )

// Modified by Ming Zeng from the original code of John Burkardt
// at https://people.math.sc.edu/Burkardt/c_src/rk4/rk4.c
/******************************************************************************/
/*
  Purpose:

    RK4VEC takes one Runge-Kutta step for a vector ODE.

  Discussion:

    It is assumed that an initial value problem, of the form

      du/dt = f ( t, u )
      u(t0) = u0

    is being solved.

    If the user can supply current values of t, u, a stepsize dt, and a
    function to evaluate the derivative, this function can compute the
    fourth-order Runge Kutta estimate to the solution at time t+dt.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T0, the current time.

    Input, int M, the spatial dimension.

    Input, double U0[M], the solution estimate at the current time.

    Input, double DT, the time step.

    Input, double *F ( double T, int M, double U[] ), a function which evaluates
    the derivative, or right hand side of the problem.

    Output, double RK4VEC[M], the fourth-order Runge-Kutta solution estimate
    at time T0+DT.
*/
{
  double *f0;
  double *f1;
  double *f2;
  double *f3;
  int i;
  double t1;
  double t2;
  double t3;
  double *u1;
  double *u2;
  double *u3;
/*
  Get four sample values of the derivative.
*/
  f0 = f ( t0, u0 );

  t1 = t0 + dt / 2.0;
  u1 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    u1[i] = u0[i] + dt * f0[i] / 2.0;
  }
  f1 = f ( t1, u1 );

  t2 = t0 + dt / 2.0;
  u2 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    u2[i] = u0[i] + dt * f1[i] / 2.0;
  }
  f2 = f ( t2, u2 );

  t3 = t0 + dt;
  u3 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
     u3[i] = u0[i] + dt * f2[i];
  }
  f3 = f ( t3, u3 );
/*
  Combine them to estimate the solution.
*/
  for ( i = 0; i < m; i++ )
  {
     u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
  }
/*
  Free memory.
*/
  free ( f0 );
  free ( f1 );
  free ( f2 );
  free ( f3 );
  free ( u1 );
  free ( u2 );
  free ( u3 );

  return 0;
}
