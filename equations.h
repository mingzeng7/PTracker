//define the forces
#ifndef __FORCE_H
#define __FORCE_H

#include "includes.h"

void f_ext ( const double * x, const double * dot_x, double * f );
/* External force.
   Input:
   x: position
      x[0] is zeta = z - beta_w*t
      x[1] is x
      x[2] is y
   dot_x: time derivative of x
      dot_x[0] is beta_z - beta_w
      dot_x[1] is beta_x
      dot_x[2] is beta_y
   Output:
   f: Pointer to an array of two doubles, containing the force.
      f[0] is f^ext_z
      f[1] is f^ext_x
      f[2] is f^ext_y
*/

void f_rad ( const double * p, const double * dp_over_dt, double * f );
/* Radiation reaction force. For 1st radiation term turned on (if_RR1=1), estimate d^2 p/dt^2 by d f^ext/dt.
   Input:
   p: momentum
      p[0] is pz
      p[1] is px
   dp_over_dt: time derivative of p
      dp_over_dt[0] is dpz/dt
      dp_over_dt[1] is dpx/dt
   Output:
   f: Pointer to an array of two doubles, containing the estimated radiation reactions.
*/

double * dy_over_dt ( double t, const double * y );
/* Create and set the array of dy_over_dt. Need to be released manually after using it. */
#endif /* #ifndef __FORCE_H */
