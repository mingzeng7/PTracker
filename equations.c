#include "equations.h"

void f_ext ( const double * x, double * f )
/* External force.
   Input:
   x: position
      x[0] is xi = z - beta_w*t
      x[1] is x
   Output:
   f: Pointer to an array of two doubles, containing the force.
      f[0] is f^ext_z
      f[1] is f^ext_x
*/
{
  f[0] = -.5*x[0];
  f[1] = -.5*x[1];
  return;
}

void f_rad ( const double * p, const double * dp_over_dt, double * f )
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
{
  const double gamma = sqrt(1.+Square(p[0])+Square(p[1]));
  const double gamma_dot = (p[0]*dp_over_dt[0]+p[1]*dp_over_dt[1])/gamma;
  if(if_RR1)
  {
    // 1st term of radiation reactions
    // d^2 p/dt^2 is replaced by d f^ext/dt.
    // We assert f^ext_x = -x/2 and f^ext_z = -zeta/2 here.
    // If f^ext is changed, remember to change the following 2 lines also.
    f[0] = re_times_2_over_3 * (gamma_dot*dp_over_dt[0] + .5*(beta_w*gamma - p[0]));
    f[1] = re_times_2_over_3 * (gamma_dot*dp_over_dt[1] - .5*p[1]);
  }
  else
  {
    f[0] = 0.;
    f[1] = 0.;
  }

  if(if_RR2)
  {
    // 2nd term of radiation reactions
    const double common_term = re_times_2_over_3 * gamma * (Square(gamma_dot) - Square(dp_over_dt[0]) - Square(dp_over_dt[1]));
    f[0] += common_term * p[0];
    f[1] += common_term * p[1];
  }
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
   Output:
      Pointer to an array of doubles, containing the results of dy_over_dt. Programmer should release this memory allocation after use.
 */
{
  double * out_buffer;
  const double gamma = sqrt(1.+Square(y[2])+Square(y[3]));
  const double RelTol = 1.e-10;//Relative Tolarance for modifying dp/dt using the RR force
  const double extreme_small = 1.e-24;//Set a extreme small number

  out_buffer = ( double * ) malloc ( 4 * sizeof * out_buffer);//programmer should release this pointer outside this function
  out_buffer[0] = y[2]/gamma - beta_w; // = beta_z
  out_buffer[1] = y[3]/gamma; // = beta_x
  // set pz dot and px dot equal to external forces first.
  f_ext(y, out_buffer+2);

  if(if_RR)
  {
    int i;
    const int max_cycles = 3;
    double f_ext[2]={out_buffer[2],out_buffer[3]};
    double tmp_f_RR[2];
    double new_tmp_f_RR[2];
    f_rad(y+2,f_ext,tmp_f_RR);//set tmp rr force with the external force only
    for(i=0;i<max_cycles;i++)
    {
      //printf("In dy_over_dt, cycle i = %d\n",i);
      out_buffer[2] = f_ext[0] + tmp_f_RR[0];
      out_buffer[3] = f_ext[1] + tmp_f_RR[1];
      
      f_rad(y+2,out_buffer+2,new_tmp_f_RR);//calculate new tmp rr force
      if((fabs(new_tmp_f_RR[0]-tmp_f_RR[0])<fabs(tmp_f_RR[0]*RelTol) || fabs(tmp_f_RR[0])<extreme_small)
         && (fabs(new_tmp_f_RR[1]-tmp_f_RR[1])<fabs(tmp_f_RR[1]*RelTol) || fabs(tmp_f_RR[1])<extreme_small))
      //Check if the change is negligible; if not, prepare for next loop
        break;
      //tmp_f_RR[0] = 0.618*new_tmp_f_RR[0]+0.382*tmp_f_RR[0];
      //tmp_f_RR[1] = 0.618*new_tmp_f_RR[1]+0.382*tmp_f_RR[0];
      tmp_f_RR[0] = new_tmp_f_RR[0];
      tmp_f_RR[1] = new_tmp_f_RR[1];
    }
    if(i>=max_cycles) printf("Warrning: max_cycles is reached in dy_over_dt!\n");
  }
  return out_buffer;
}
