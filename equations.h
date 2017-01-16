//define the forces
#ifndef __FORCE_H
#define __FORCE_H

#include "includes.h"

void f_rad2 ( const double * p, const double * dp_over_dt, double * f );
/* Radiation reaction force of the second term.
   u is an double array with 4 elements
   p[0] is the longitudinal momentum p1
   p[1] is the tramsverse momentum p2
   dp_over_dt[0] is dp1/dt
   dp_over_dt[1] is dp2/dt */

double * dy_over_dt ( double t, const double * y );
/* Create and set the array of dy_over_dt. Need to be released manually after using it. */
#endif /* #ifndef __FORCE_H */
