/* Constants definition header file*/
#ifndef __CONSTANT_H
#define __CONSTANT_H

#define Square(x) ((x)*(x))
//default information for the input file
extern const int	max_file_name_length;
extern const char*	default_input_file;

//default information for the output file
extern const char*      default_output_file;

//default values for the input
extern const unsigned short int default_allocation_method;
extern const unsigned long int	max_particle_number;
extern const unsigned short int	default_particle_number;

//radomness
extern const int 		default_time_seed;
extern const unsigned int	default_seed;

//physical parameters
extern const int	default_if_RR; //default bool for whether calculate radiation reaction
extern const int	default_if_RR1; //default bool for whether calculate radiation reaction term 1
extern const int	default_if_RR2; //default bool for whether calculate radiation reaction term 2


extern const double	re_SI;
extern const double	default_kp_SI;
extern const double	default_gamma_w;
extern const double	default_f_z0;
extern const double	default_kappa_square;
extern const double	default_lambda;
extern const double	default_t0;
extern const double	default_t1;
extern const double	default_sampling_factor;//

extern const double h_bar_SI;
extern const double elec_anomaly;
extern const double elec_mc_si;
extern const double s0;
#endif /* ifndef __CONSTANT_H */
