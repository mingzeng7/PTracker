/* ----------------------------------------------------------------------------
   PTracker - by Ming Zeng, Email: chippen8888@gmail.com
   ----------------------------------------------------------------------------
*/

#include "includes.h"
int if_RR; // if_RR == 0 will turn off all the radiation reactions
int if_RR1; // if_RR1 == 0 will turn off the 1st term of radiation reactions
int if_RR2; // if_RR2 == 0 will turn off the 2nd term of radiation reactions
double re_times_2_over_3; // Classical electron radius normalized to k_p^{-1}, times 2/3
double beta_w;

int main(int argc, char **argv)
{
  const char * input_file;
  const char * output_file;
  config_t cfg;
  config_setting_t *setting;
  unsigned int particle_number;

  if(2>argc)
  {
    input_file = default_input_file;
  }
  else
    input_file = argv[1];
  printf("The input file name is: %s\n",(const char*)input_file);
  config_init(&cfg);

  /* Read the file. If there is an error, report it and exit. */
  if(! config_read_file(&cfg, (const char*)input_file))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    return(EXIT_FAILURE);
  }

  //Reading and set if_RR
  if(config_lookup_bool(&cfg, "if_RR", &if_RR))
  {
    printf("if_RR =  %d\n", if_RR);
  }
  else
  {
    if_RR = default_if_RR;
    printf("if_RR =  %d using the default value!\n", if_RR);
  }

  //Reading and set if_RR1
  if(config_lookup_bool(&cfg, "if_RR1", &if_RR1))
  {
    printf("if_RR1 =  %d\n", if_RR1);
  }
  else
  {
    if_RR1 = default_if_RR1;
    printf("if_RR1 =  %d using the default value!\n", if_RR1);
  }

  //Reading and set if_RR2
  if(config_lookup_bool(&cfg, "if_RR2", &if_RR2))
  {
    printf("if_RR2 =  %d\n", if_RR2);
  }
  else
  {
    if_RR2 = default_if_RR2;
    printf("if_RR2 =  %d using the default value!\n", if_RR2);
  }

  //Reading kp and set re
  double kp_SI;
  if(config_lookup_float(&cfg, "wake.kp_SI", &kp_SI))
  {
    printf("wake.kp_SI =  %f\n", kp_SI);
  }
  else
  {
    kp_SI = default_kp_SI;
    printf("wake.kp_SI = %f using the default value!\n", kp_SI);
  }
  re_times_2_over_3 = kp_SI * re_SI * 2.0/3.0;

  //Reading gamma_w and set beta_w
  double gamma_w;
  if(config_lookup_float(&cfg, "wake.gamma_w", &gamma_w))
  {
    printf("wake.gamma_w = %f\n", gamma_w);
  }
  else
  {
    gamma_w = default_gamma_w;
    printf("wake.gamma_w = %f using the default value!\n", gamma_w);
  }
  beta_w = sqrt(1.-1./Square(gamma_w));
  printf("beta_w = %.*e\n", 10, beta_w);

  //randomness
  int time_seed;
  int seed;
  if(config_lookup_bool(&cfg, "random.time_seed", &time_seed))
  {
    printf("random.time_seed = %d\n", time_seed);
  }
  else
  {
    time_seed = default_time_seed;
    printf("random.time_seed = %d using the default value!\n", time_seed);
  }
  if(time_seed)
  {
    srand(time(NULL));
  }
  else
  {
    if(config_lookup_int(&cfg, "random.seed", &seed))
    {
      printf("random.seed = %d\n", seed);
    }
    else
    {
      seed = default_seed;
      printf("random.seed = %d using the default value!\n", seed);
    }
    srand(seed);
  }
  //printf("random numbers %d, %d, %d, %d\n",rand(),rand(),rand(),rand());
  /* Get the output file name. */
  if(config_lookup_string(&cfg, "output_file", &output_file))
  {
    printf("Output file name: %s\n\n", output_file);
  }
  else
  {
    output_file = default_output_file;
    printf("No 'output_file' specified. Use %s as a default!\n",output_file);
  }

  double t0;
  double t1;
  if(config_lookup_float(&cfg, "t0",&t0))
    fprintf(stdout,"t0 = %f\n", t0);
  else
  {
    t0 = default_t0;
    fprintf(stdout,"t0 = %f by default!\n", t0);
  }
  if(config_lookup_float(&cfg, "t1",&t1))
    fprintf(stdout,"t1 = %f\n", t1);
  else
  {
    t1 = default_t1;
    fprintf(stdout,"t1 = %f by default!\n", t1);
  }
  double t_duration = t1 - t0;

  double sampling_factor;
  if(config_lookup_float(&cfg, "sampling_factor",&sampling_factor))
    fprintf(stdout, "sampling_factor = %f\n", sampling_factor);
  else
  {
    sampling_factor = default_sampling_factor;
    fprintf(stdout, "sampling_factor = %f by default!\n", sampling_factor);
  }

  unsigned int allocation_method=0;
  if(config_lookup_int(&cfg, "particle.allocation_method", &allocation_method))
    fprintf(stdout,"allocation_method = %d\n", allocation_method);
  else
    fprintf(stdout,"allocation_method = %d by default!\n", allocation_method);

  //start h5 creationg
  hid_t ofile_h5id;
  hid_t dset_h5id;
  hid_t space_h5id;
  herr_t h5status;
  ofile_h5id = H5Fcreate (output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  switch(allocation_method){
    default: //single particle mode
      setting = config_lookup(&cfg, "particle.single_particle");
      if(NULL == setting)
      {
        fprintf(stderr, "No single_particle set assinged in the single particle mode!\n");
        config_destroy(&cfg);
        return(EXIT_FAILURE);
      }
      else
      {
        particle_number = config_setting_length(setting);
        int i;
        char i_str[20]; // buffer for converting i to string
        for(i = 0; i < particle_number; ++i)
        {
          printf("Single particle 1:\n");
          config_setting_t * single_particle = config_setting_get_elem(setting, i);
          double A, B, gamma0;
          config_setting_lookup_float(single_particle, "A", &A);
          config_setting_lookup_float(single_particle, "B", &B);
          config_setting_lookup_float(single_particle, "gamma0", &gamma0);
          double r_omega_beta = sqrt(gamma0*2.); // 1 over Betatron oscillation frequency
          printf("Start to track particle #%d: A = %f, B = %f, gamma0 = %f\n", i,A,B, gamma0);
          double dt = r_omega_beta/sampling_factor;
          fprintf(stdout,"dt = %f\n", dt);
          int track_arry_lenth = ((int) (t_duration / dt))+1;
          fprintf(stdout,"track_arry_lenth =  %d\n", track_arry_lenth);
          double * buffer;
          int column_length = 7; // The columns are t, zeta, x, y, pz, px, py
          hsize_t * dims;
          dims=(hsize_t *)malloc(2*sizeof *dims);
          dims[0]=track_arry_lenth;
          dims[1]=column_length;

          buffer = (double*)malloc(dims[1]*dims[0]*sizeof *buffer);
          buffer[0] = t0;

          // Get intial phase of the oscillation
          double phase0;
          if(config_setting_lookup_float(single_particle, "phase0", &phase0))
          {
            fprintf(stdout,"Particle # %d phase0 = %f\n", i, phase0);
            phase0 *= M_PI/180; // Transform degree to radian
          }
          else
          {
            fprintf(stderr,"Particle # %d phase0 not given! Set to default 0.\n", i);
            phase0 = 0.;
          }

          double x_1 = A/sqrt(sqrt(gamma0)); // Amplitude of x oscillation
          double y_1 = B/sqrt(sqrt(gamma0)); // Amplitude of y oscillation
          double zeta0;
          if(config_setting_lookup_float(single_particle, "zeta0", &zeta0))
            fprintf(stdout,"Particle # %d zeta0 = %f\n", i, zeta0);
          else
            fprintf(stderr,"Particle # %d zeta0 not given!\n", i);
          // Set buffer[1] to be initial value of zeta = zeta0 with phase shift
          buffer[1] = zeta0 - (Square(x_1)+Square(y_1))/r_omega_beta/8*sin(phase0*2);

          // Set buffer[2] to be initial value of x
          buffer[2] = x_1*sin(phase0);
          // Set buffer[3] to be initial value of y
          buffer[3] = y_1*cos(phase0);
          // Set buffer[4] to be initial value of pz with phase shift
          buffer[4] = sqrt(Square(gamma0)-1.0-gamma0*(Square(x_1)+Square(y_1))*(4+cos(phase0*2))/16);
          // Set buffer[5] to be initial value of px with phase shift
          buffer[5] = gamma0*x_1*cos(phase0)/r_omega_beta;
          // Set buffer[6] to be initial value of py with phase shift
          buffer[6] = gamma0*y_1*sin(phase0)/r_omega_beta;

          // j is the index for rows
          int j;
          for(j = 1; j < track_arry_lenth; j++)
          {
            int pre_row_start_index = (j-1)*column_length;
            int this_row_start_index = j*column_length;
            buffer[this_row_start_index] = buffer[pre_row_start_index] + dt;
            rk4vec ( buffer[pre_row_start_index], column_length-1, buffer+pre_row_start_index+1, buffer+this_row_start_index+1, dt, dy_over_dt );
          }
          sprintf(i_str, "%d", i);
          h5status = H5LTmake_dataset_double (ofile_h5id, i_str, 2, dims ,(const double*)buffer);
          free(buffer);
          free(dims);
        }
      }
  }

  config_destroy(&cfg);
  h5status = H5Fclose (ofile_h5id);
  printf("h5status = %d\n", h5status);
  return(EXIT_SUCCESS);
}

/* eof */
