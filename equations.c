#include "equations.h"

void f_rad2 ( const double * p, const double * dp_over_dt, double * f )
/* Radiation reaction force of the second term.
   u is an double array with 4 elements
   p[0] is the longitudinal momentum p1
   p[1] is the tramsverse momentum p2
   dp_over_dt[2] is dp1/dt
   dp_over_dt[3] is dp2/dt */
{
  const double gamma = sqrt(1.+Square(p[0])+Square(p[1]));
  const double gamma_dot = p[0]/gamma*dp_over_dt[0]+p[1]/gamma*dp_over_dt[1];
  const double common_term = re_times_2_over_3 * gamma * (Square(gamma_dot) - Square(dp_over_dt[0]) - Square(dp_over_dt[1]));
  f[0] = common_term * p[0];
  f[1] = common_term * p[1];
  return;
}

double * dy_over_dt ( double t, const double * y )
/* Create and set the array of dy_overe_dt. Programmer should release this pointer after use.
 * y[0] is z-beta_w * t, where z is the longitudinal position, and beta_w is the speed of the wake.
 * y[1] is x.
 * y[2] is pz.
 * y[3] ps px.
 */
{
  double * out_buffer;
  const double gamma = sqrt(1.+Square(y[2])+Square(y[3]));
  const double RelTol = 1.e-6;//Relative Tolarance for modifying dp/dt using the RR force
  const double extreme_small = 1.e-24;//Set a extreme small number

  out_buffer = ( double * ) malloc ( 4 * sizeof * out_buffer);//programmer should release this pointer outside this function
  out_buffer[0] = y[2]/gamma - beta_w;
  out_buffer[1] = y[3]/gamma;
  // set p1 dot and p2 dot with external force first.
  out_buffer[2] = y[0]/(-2.);
  out_buffer[3] = y[1]/(-2.);
  //out_buffer[3] = y[1]/(-4.)*(1.+y[2]/gamma); // The actual magnetic force is a bit smaller than electric force

  int i;
  if(if_RR)
  {
    const int max_cycles = 10;
    double f_ext[2]={out_buffer[2],out_buffer[3]};
    double tmp_f_RR[2];
    double new_tmp_f_RR[2];
    f_rad2(y+2,f_ext,tmp_f_RR);//set tmp rr force with the external force only
    for(i=0;i<max_cycles;i++)
    {
      //printf("In dy_over_dt, cycle i = %d\n",i);
      out_buffer[2] = f_ext[0] + tmp_f_RR[0];
      out_buffer[3] = f_ext[1] + tmp_f_RR[1];
      
      f_rad2(y+2,out_buffer+2,new_tmp_f_RR);//calculate new tmp rr force
      if((fabs(new_tmp_f_RR[0]-tmp_f_RR[0])<fabs(tmp_f_RR[0]*RelTol) || fabs(tmp_f_RR[0])<extreme_small)
         && (fabs(new_tmp_f_RR[1]-tmp_f_RR[1])<fabs(tmp_f_RR[1]*RelTol) || fabs(tmp_f_RR[1])<extreme_small))
      //Check if the change is large enough; if not, break
        break;
      //printf("f_RR[0] is changed from %.*e to %.*e; relative change is %.*e\n",20,tmp_f_RR[0],20,new_tmp_f_RR[0],10,(new_tmp_f_RR[0]-tmp_f_RR[0])/tmp_f_RR[0]);
      //printf("f_RR[1] is changed from %.*e to %.*e; relative change is %.*e\n",20,tmp_f_RR[1],20,new_tmp_f_RR[1],10,(new_tmp_f_RR[1]-tmp_f_RR[1])/tmp_f_RR[1]);
      //tmp_f_RR[0] = 0.618*new_tmp_f_RR[0]+0.382*tmp_f_RR[0];
      tmp_f_RR[0] = new_tmp_f_RR[0];
      //tmp_f_RR[1] = 0.618*new_tmp_f_RR[1]+0.382*tmp_f_RR[0];
      tmp_f_RR[1] = new_tmp_f_RR[1];
    }
    if(i>=max_cycles) printf("Warrning: i>=max_cycles in dy_over_dt!\n");
  }
  return out_buffer;
}
