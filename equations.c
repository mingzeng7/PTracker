#include "equations.h"

void f_rad2 ( const double * p, const double * dp_over_dt, double * f )
/* Radiation reaction force of the second term.
   Input:
   p: momentum
      p[0] is pz
      p[1] is px
   dp_over_dt: time derivative of p
      dp_over_dt[0] is dpz/dt
      dp_over_dt[1] is dpx/dt
   Output:
   f: Pointer to an array of two doubles, containing the estimated 2nd term of the radiation reaction
*/
{
  const double gamma = sqrt(1.+Square(p[0])+Square(p[1]));
  const double gamma_dot = p[0]/gamma*dp_over_dt[0]+p[1]/gamma*dp_over_dt[1];
  const double common_term = re_times_2_over_3 * gamma * (Square(gamma_dot) - Square(dp_over_dt[0]) - Square(dp_over_dt[1]));
  f[0] = common_term * p[0];
  f[1] = common_term * p[1];
  return;
}

double * dy_over_dt ( double t, const double * y )
/* Create and set the array of dy_over_dt.
   Input:
   t: time
   y: y[0] is z-beta_w * t, where z is the longitudinal position, and beta_w is the speed of the wake.
      y[1] is x.
      y[2] is pz.
      y[3] ps px.
      y[4] is dpz/dt.
      y[5] ps dpx/dt.
   Output:
      Pointer to an array of doubles, containing the results of dy_over_dt. Programmer should release this memory allocation after use.
 */
{
  double * out_buffer;
  const double gamma = sqrt(1.+Square(y[2])+Square(y[3]));
  const double beta_z = y[2]/gamma;

  out_buffer = ( double * ) malloc ( 6 * sizeof * out_buffer);//programmer should release this pointer outside this function
  out_buffer[0] = beta_z - beta_w;
  out_buffer[1] = y[3]/gamma; // out_buffer[1] is beta_x

  if(if_RR)
  {  
    const double dgamma_over_dt = beta_z*y[4] + out_buffer[1]*y[5];
    out_buffer[2] = y[4];
    out_buffer[3] = y[5];
    const double tmp = Square(dgamma_over_dt)-Square(y[4])-Square(y[5]);
    out_buffer[4] = ((y[4]+y[0]*.5)/re_times_2_over_3 - dgamma_over_dt*y[4])/gamma - tmp*y[2];
    out_buffer[5] = ((y[5]+y[1]*.5)/re_times_2_over_3 - dgamma_over_dt*y[5])/gamma - tmp*y[3];
    fprintf(stdout,"%e %e %e %e %e %e\n", y[0], y[1], y[2], y[3], y[4], y[5]);
  }
  else
  {
    // Do not have to calculate out_buffer[4] and out_buffer[5] if the radiation reaction is turned off
    out_buffer[2] = -0.5*y[0];
    out_buffer[3] = -0.5*y[1];
  }
  return out_buffer;
}
