#include "equations.h"

void f_ext ( const double * x, const double * dot_x, double * f )
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
   f: Pointer to an array of 3 doubles, containing the force.
      f[0] is f^ext_z = f_z0 - lambda*zeta + kappa_square*lambda*(x*beta_x + y*beta_y);
      f[1] is f^ext_x = -kappa_square*(1.-lambda+lambda*beta_z)*x
      f[2] is f^ext_y = -kappa_square*(1.-lambda+lambda*beta_z)*y
*/
{
  const double common_term = -kappa_square*(1.+lambda*(beta_w-1.+dot_x[0]));
  f[0] = f_z0 - lambda*x[0] + kappa_square*lambda*(x[1]*dot_x[1]+x[2]*dot_x[2]);
  f[1] = x[1]*common_term;
  f[2] = x[2]*common_term;
  return;
}

void f_rad ( const double * p, const double * dp_over_dt, double * f )
/* Radiation reaction force. For 1st radiation term turned on (if_RR1=1), estimate d^2 p/dt^2 by d f^ext/dt.
   Input:
   p: momentum
      p[0] is pz
      p[1] is px
      p[2] is py
   dp_over_dt: time derivative of p
      dp_over_dt[0] is dpz/dt
      dp_over_dt[1] is dpx/dt
      dp_over_dt[2] is dpy/dt
   Output:
   f: Pointer to an array of 3 doubles, containing the estimated radiation reactions.
*/
{
  const double gamma = sqrt(1.+Square(p[0])+Square(p[1])+Square(p[2]));
  const double gamma_dot = (p[0]*dp_over_dt[0]+p[1]*dp_over_dt[1]+p[2]*dp_over_dt[2])/gamma;
  if(if_RR1)
  {
    // 1st term of radiation reactions
    // d^2 p/dt^2 is replaced by d f^ext/dt.
    // We assert f^ext_x = -x/2 and f^ext_z = -zeta/2 here.
    // If f^ext changes, remember to change the following 3 lines also.
    f[0] = re_times_2_over_3 * (gamma_dot*dp_over_dt[0] + .5*(beta_w*gamma - p[0]));
    f[1] = re_times_2_over_3 * (gamma_dot*dp_over_dt[1] - .5*p[1]);
    f[2] = re_times_2_over_3 * (gamma_dot*dp_over_dt[2] - .5*p[2]);
  }
  else
  {
    f[0] = 0.;
    f[1] = 0.;
    f[2] = 0.;
  }

  if(if_RR2)
  {
    // 2nd term of radiation reactions
    const double common_term = re_times_2_over_3 * gamma * (Square(gamma_dot) - Square(dp_over_dt[0]) - Square(dp_over_dt[1]) - Square(dp_over_dt[2]));
    f[0] += common_term * p[0];
    f[1] += common_term * p[1];
    f[2] += common_term * p[2];
  }
  return;
}

double * dy_over_dt ( double t, const double * y )
/* Create and set the array of dy_over_dt.
   Input:
   t: time
   y: y[0] is zeta = z-beta_w * t, where z is the longitudinal position, and beta_w is the speed of the wake.
      y[1] is x.
      y[2] is y.
      y[3] is pz.
      y[4] ps px.
      y[5] ps py.
   Output:
      Pointer to an array of doubles, containing the results of dy_over_dt. Programmer should release this memory allocation after use.
 */
{
  double * out_buffer;
  const double gamma = sqrt(1.+Square(y[3])+Square(y[4])+Square(y[5]));
  const double RelTol = 1.e-5;//Relative Tolarance for modifying dp/dt using the RR force
  const double extreme_small = 1.e-10;//Set an extremely small number

  out_buffer = ( double * ) malloc ( 6 * sizeof * out_buffer);//programmer should release this pointer outside this function
  out_buffer[0] = y[3]/gamma - beta_w; // = beta_z - beta_w = dot_zeta
  out_buffer[1] = y[4]/gamma; // = beta_x
  out_buffer[2] = y[5]/gamma; // = beta_y
  // set pz dot and px dot equal to external forces first.
  f_ext(y, out_buffer, out_buffer+3);

  if(if_RR)
  {
    int i;
    const int max_cycles = 3;
    double f_ext_save[3]={out_buffer[3],out_buffer[4],out_buffer[5]};
    double tmp_f_RR[3];
    double new_tmp_f_RR[3];
    double diff_tmp_f_RR[3];
    f_rad(y+3,f_ext_save,tmp_f_RR);//set tmp rr force with the external force only
    for(i=0;i<max_cycles;i++)
    {
      //printf("In dy_over_dt, cycle i = %d\n",i);
      out_buffer[3] = f_ext_save[0] + tmp_f_RR[0];
      out_buffer[4] = f_ext_save[1] + tmp_f_RR[1];
      out_buffer[5] = f_ext_save[2] + tmp_f_RR[2];
      
      f_rad(y+3,out_buffer+3,new_tmp_f_RR);//calculate new tmp rr force
      diff_tmp_f_RR[0] = fabs(new_tmp_f_RR[0]-tmp_f_RR[0]);
      diff_tmp_f_RR[1] = fabs(new_tmp_f_RR[1]-tmp_f_RR[1]);
      diff_tmp_f_RR[2] = fabs(new_tmp_f_RR[2]-tmp_f_RR[2]);
      if((diff_tmp_f_RR[0]<fabs(tmp_f_RR[0]*RelTol) || diff_tmp_f_RR[0]<extreme_small)
         && (diff_tmp_f_RR[1]<fabs(tmp_f_RR[1]*RelTol) || diff_tmp_f_RR[1]<extreme_small)
         && (diff_tmp_f_RR[2]<fabs(tmp_f_RR[2]*RelTol) || diff_tmp_f_RR[2]<extreme_small))
      //Check if the change is negligible; if not, prepare for next loop
        break;
      //tmp_f_RR[0] = 0.618*new_tmp_f_RR[0]+0.382*tmp_f_RR[0];
      //tmp_f_RR[1] = 0.618*new_tmp_f_RR[1]+0.382*tmp_f_RR[0];
      tmp_f_RR[0] = new_tmp_f_RR[0];
      tmp_f_RR[1] = new_tmp_f_RR[1];
      tmp_f_RR[2] = new_tmp_f_RR[2];
    }
    if(i>=max_cycles) printf("Warrning: max_cycles is reached in dy_over_dt!\n");
  }
  return out_buffer;
}
