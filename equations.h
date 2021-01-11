//define the forces
#ifndef __FORCE_H
#define __FORCE_H

#include "includes.h"

void f_rad2 ( const double * p, const double * dp_over_dt, double * f );
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

double * dy_over_dt ( double t, const double * y );
/* Create and set the array of dy_over_dt. Need to be released manually after using it. */
#endif /* #ifndef __FORCE_H */
