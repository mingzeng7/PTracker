/* ----------------------------------------------------------------------------
   PTracker - by Ming Zeng, Email: chippen8888@gmail.com
   ----------------------------------------------------------------------------
*/

#include "includes.h"
#include <omp.h>  

int if_RR; // if_RR == 0 will turn off all the radiation reactions
int if_RR1; // if_RR1 == 0 will turn off the 1st term of radiation reactions
int if_RR2; // if_RR2 == 0 will turn off the 2nd term of radiation reactions
double re_times_2_over_3; // Classical electron radius normalized to k_p^{-1}, times 2/3
double beta_w; // The wake phase velocity normalized to c
double f_z0; // A constant external force in the z direction
double lambda; // The slop of longitudinal electric field, E_z = lambda * zeta
double kappa_square; // The transverse restoring parameter, E_r = kappa_square * (1-lambda) * r, B_theta = -kappa_square * lambda * r
double kappa_square_lambda; // = kappa_square * lambda
double h_bar_over_2;//normalized h_bar/2

// Gaussian random particle generation

double gaussian_rand_threadsafe(double mean, double stddev, unsigned int *seed) {
    double u, v, z;
    
    u = (rand_r(seed) + 1.0) / (RAND_MAX + 2.0);
    v = rand_r(seed) / (RAND_MAX + 1.0);
    z = sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
    
    return mean + stddev * z;
}
// Uniform random generation of particles

double uniform_rand_threadsafe(double min, double max, unsigned int *seed) {
    return min + (max - min) * (rand_r(seed) / (double)RAND_MAX);
}

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

  //Read and set if_RR
  if(config_lookup_bool(&cfg, "if_RR", &if_RR))
  {
    printf("if_RR =  %d\n", if_RR);
  }
  else
  {
    if_RR = default_if_RR;
    printf("if_RR =  %d using the default value!\n", if_RR);
  }

  //Read and set if_RR1
  if(config_lookup_bool(&cfg, "if_RR1", &if_RR1))
  {
    printf("if_RR1 =  %d\n", if_RR1);
  }
  else
  {
    if_RR1 = default_if_RR1;
    printf("if_RR1 =  %d using the default value!\n", if_RR1);
  }

  //Read and set if_RR2
  if(config_lookup_bool(&cfg, "if_RR2", &if_RR2))
  {
    printf("if_RR2 =  %d\n", if_RR2);
  }
  else
  {
    if_RR2 = default_if_RR2;
    printf("if_RR2 =  %d using the default value!\n", if_RR2);
  }

  //Read kp and set re and h_bar
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
  h_bar_over_2 = h_bar_SI * kp_SI / 2/elec_mc_si;
  printf("hbar is:%e\n", h_bar_over_2);

  //Read gamma_w and set beta_w
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

  //Read and set f_z0
  if(config_lookup_float(&cfg, "wake.f_z0", &f_z0))
  {
    printf("wake.f_z0 =  %f\n", f_z0);
  }
  else
  {
    f_z0 = default_f_z0;
    printf("wake.f_z0 = %f using the default value!\n", f_z0);
  }

  //Read and set kappa_square
  if(config_lookup_float(&cfg, "wake.kappa_square", &kappa_square))
  {
    printf("wake.kappa_square =  %f\n", kappa_square);
  }
  else
  {
    kappa_square = default_kappa_square;
    printf("wake.kappa_square = %f using the default value!\n", kappa_square);
  }

  //Read and set lambda
  if(config_lookup_float(&cfg, "wake.lambda", &lambda))
  {
    printf("wake.lambda =  %f\n", lambda);
  }
  else
  {
    lambda = default_lambda;
    printf("wake.lambda = %f using the default value!\n", lambda);
  }

  //Set kappa_square_lambda
  kappa_square_lambda = kappa_square * lambda;

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
    seed = time(NULL); 
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
  
  // Set the number of OpenMP threads
  int num_threads = 0;
  if(config_lookup_int(&cfg, "parallel.num_threads", &num_threads))
  {
    printf("Using %d OpenMP threads\n", num_threads);
    omp_set_num_threads(num_threads);
  }
  else
  {
    printf("Using default number of OpenMP threads\n");
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
  //The switch function executes different branch code blocks based on the value of 'allocation.method'
  switch(allocation_method){
    case 0://single particle mode
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
          double x_1; // Amplitude of x oscillation. Setting x_1 overrides A
          double y_1; // Amplitude of y oscillation. Setting y_1 overrides B
          
          config_setting_lookup_float(single_particle, "gamma0", &gamma0);
          if(config_setting_lookup_float(single_particle, "x_1", &x_1))
          {
            fprintf(stdout,"Particle # %d x_1 = %f\n", i, x_1);
            A = x_1 * sqrt(sqrt(gamma0)); // "A" may be not necessary. This line may be deleted.
          }
          else
          {
            config_setting_lookup_float(single_particle, "A", &A);
            x_1 = A/sqrt(sqrt(gamma0));
          }
          if(config_setting_lookup_float(single_particle, "y_1", &y_1))
          {
            fprintf(stdout,"Particle # %d y_1 = %f\n", i, y_1);
            B = y_1 * sqrt(sqrt(gamma0)); // "B" may be not necessary. This line may be deleted.
          }
          else
          {
            config_setting_lookup_float(single_particle, "B", &B);
            y_1 = B/sqrt(sqrt(gamma0));
          }
          
          double r_omega_beta = sqrt(gamma0/kappa_square); // 1 over Betatron oscillation frequency
          printf("Start to track particle #%d: x_1 = %f, y_1 = %f, gamma0 = %f\n", i, x_1, y_1, gamma0);
          double dt = r_omega_beta/sampling_factor;
          fprintf(stdout,"dt = %f\n", dt);
          int track_arry_lenth = ((int) (t_duration / dt))+1;
          fprintf(stdout,"track_arry_lenth =  %d\n", track_arry_lenth);
          double * buffer;
          int column_length = 10; // The columns are t, zeta, x, y, pz, px, py sz,sx,sy,
          hsize_t * dims;
          dims=(hsize_t *)malloc(2*sizeof *dims);
          dims[0]=track_arry_lenth;
          dims[1]=column_length;

          buffer = (double*)malloc(dims[1]*dims[0]*sizeof *buffer);
          buffer[0] = t0;

          // Get intial phase of the oscillation
          double phase0_x, phase0_y;
          if(config_setting_lookup_float(single_particle, "phase0_x", &phase0_x))
          {
            fprintf(stdout,"Particle # %d phase0_x = %f\n degree", i, phase0_x);
            phase0_x *= M_PI/180; // Transform degree to radian
          }
          else
          {
            fprintf(stderr,"Particle # %d phase0_x not given! Set to default 0.\n", i);
            phase0_x = 0.;
          }
          if(config_setting_lookup_float(single_particle, "phase0_y", &phase0_y))
          {
            fprintf(stdout,"Particle # %d phase0_y = %f\n degree", i, phase0_y);
            phase0_y *= M_PI/180; // Transform degree to radian
          }
          else
          {
            fprintf(stderr,"Particle # %d phase0_y not given! Set to default 0.\n", i);
            phase0_y = 0.;
          }
          //get initial s phase
          double theta0, phi0;
          if(config_setting_lookup_float(single_particle, "theta0", &theta0))
          {
            fprintf(stdout,"Particle # %d theta0 = %f\n degree", i, theta0);
            theta0 *= M_PI/180; // Transform degree to radian
          }
          else
          {
            fprintf(stderr,"Particle # %d theta0 not given! Set to default 0.\n", i);
            theta0 = 0.;
          }
          if(config_setting_lookup_float(single_particle, "phi0", &phi0))
          {
            fprintf(stdout,"Particle # %d phi0 = %f\n degree", i, phi0);
            phi0 *= M_PI/180; // Transform degree to radian
          }
          else
          {
            fprintf(stderr,"Particle # %d phi0 not given! Set to default 0.\n", i);
            phi0 = 0.;
          }
          // The intial betatron oscillation follows
          // zeta = zeta0 - (x * beta_x + y * beta_y)/4,
          // x = x_1 * cos (omega_beta * t + phase0_x)
          // y = y_1 * cos (omega_beta * t + phase0_y)
          // beta_x = -omega_beta * x_1 * sin (omega_beta * t + phase0_x)
          // beta_y = -omega_beta * y_1 * sin (omega_beta * t + phase0_y)
          // px = -gamma0 * omega_beta * x_1 * sin (omega_beta * t + phase0_x)
          // py = -gamma0 * omega_beta * y_1 * sin (omega_beta * t + phase0_y)

          double zeta0;
          if(config_setting_lookup_float(single_particle, "zeta0", &zeta0))
            fprintf(stdout,"Particle # %d zeta0 = %f\n", i, zeta0);
          else
            fprintf(stderr,"Particle # %d zeta0 not given!\n", i);
          // Set buffer[1] to be initial value of zeta = zeta0 - 0.25*(x*beta_x + y*beta_y)
          buffer[1] = zeta0 + (Square(x_1)*sin(phase0_x*2)+Square(y_1)*sin(phase0_y*2))/r_omega_beta/8;

          // Set buffer[2] to be initial value of x
          buffer[2] = x_1*cos(phase0_x);
          // Set buffer[3] to be initial value of y
          buffer[3] = y_1*cos(phase0_y);
          // Set buffer[5] to be initial value of px
          buffer[5] = -gamma0*x_1*sin(phase0_x)/r_omega_beta;
          // Set buffer[6] to be initial value of py
          buffer[6] = -gamma0*y_1*sin(phase0_y)/r_omega_beta;
          // Set buffer[4] to be initial value of pz
          // For precise initial value we need beta_z0
          double beta_z0 = 1. - 0.5*(1./Square(gamma0) + 0.5*kappa_square/gamma0*(Square(x_1)+Square(y_1)));
          // gamma = gamma0 + 0.5*(kappa_square*(1.-lambda) - lambda*beta_z0*0.25) * ((Square(x_1) + Square(y_1))/2 - Square(buffer[2])-Square(buffer[3]))
          buffer[4] = sqrt(Square(gamma0 + 0.5*(kappa_square*(1.-lambda) - lambda*beta_z0*0.25) * ((Square(x_1) + Square(y_1))/2 - Square(buffer[2])-Square(buffer[3])))-1.-Square(buffer[5])-Square(buffer[6]));
          //Set buffer[7] to be initial value of sz
          buffer[7] = s0 * cos(theta0);
          //set buffer[8] to be initial value of sx
          buffer[8] = s0 * sin(theta0)*cos(phi0);
          //set buffer[9] to be initial value of sy
          buffer[9] = s0 * sin(theta0) * sin(phi0);
          // j is the index for rows
          int j;
          for(j = 1; j < track_arry_lenth; j++)
          {
            int pre_row_start_index = (j-1)*column_length;
            int this_row_start_index = j*column_length;
            buffer[this_row_start_index] = buffer[pre_row_start_index] + dt;
            int opt = rk4vec ( buffer[pre_row_start_index], column_length-1, buffer+pre_row_start_index+1, buffer+this_row_start_index+1, dt, dy_over_dt );
          }
          sprintf(i_str, "%d", i);
          h5status = H5LTmake_dataset_double (ofile_h5id, i_str, 2, dims ,(const double*)buffer);
          free(buffer);
          free(dims);
        }
      }
      break;
    case 1: // Multi-particle bunch mode
      {
          printf("Multi-particle bunch mode selected (batch processing to save memory).\n");
          
          // Read bunch parameters
          config_setting_t *bunch_setting = config_lookup(&cfg, "particle.bunch");
          if(NULL == bunch_setting)
          {
              fprintf(stderr, "No bunch parameters specified in multi-particle mode!\n");
              config_destroy(&cfg);
              return(EXIT_FAILURE);
          }
          
          // Read bunch size
          int bunch_size;
          if(!config_setting_lookup_int(bunch_setting, "number", &bunch_size))
          {
              fprintf(stderr, "No bunch size specified! Using default 100.\n");
              bunch_size = 100;
          }
          printf("Bunch size: %d\n", bunch_size);
          
          // Read distribution type
          const char* distribution_type = "gaussian";
          config_setting_lookup_string(bunch_setting, "distribution", &distribution_type);
          printf("Distribution type: %s\n", distribution_type);
          
          // Read storage sampling factor for bunch mode
          double storage_sampling_factor = 1.0; // each time step is stored by default
          if(config_setting_lookup_float(bunch_setting, "storage_sampling_factor", &storage_sampling_factor))
          {
              printf("Storage sampling factor = %f\n", storage_sampling_factor);
          }
          else
          {
              storage_sampling_factor = 1.0;
              printf("Using default storage sampling factor = 1.0 (store every time step)\n");
          }
          
          if(storage_sampling_factor < 1.0) {
              storage_sampling_factor = 1.0;
              printf("Warning: storage_sampling_factor cannot be less than 1.0, using 1.0\n");
          }
          
          // Read bunch parameters for each dimension
          double zeta0_mean, zeta0_std;
          double x_mean, x_std, y_mean, y_std;
          double px_mean, px_std, py_mean, py_std;
          double gamma0_mean, gamma0_std;
          double theta0_mean, theta0_std, phi0_mean, phi0_std;
          
          // Longitudinal position (zeta)
          if(!config_setting_lookup_float(bunch_setting, "zeta0_mean", &zeta0_mean))
              zeta0_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "zeta0_std", &zeta0_std))
              zeta0_std = 0.1;
          
          // Transverse positions
          if(!config_setting_lookup_float(bunch_setting, "x_mean", &x_mean))
              x_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "x_std", &x_std))
              x_std = 1.0e-6;
          
          if(!config_setting_lookup_float(bunch_setting, "y_mean", &y_mean))
              y_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "y_std", &y_std))
              y_std = 1.0e-6;
          
          // Transverse momenta
          if(!config_setting_lookup_float(bunch_setting, "px_mean", &px_mean))
              px_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "px_std", &px_std))
              px_std = 0.0;
          
          if(!config_setting_lookup_float(bunch_setting, "py_mean", &py_mean))
              py_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "py_std", &py_std))
              py_std = 0.0;
          
          // Energy (gamma0)
          if(!config_setting_lookup_float(bunch_setting, "gamma0_mean", &gamma0_mean))
              gamma0_mean = 1.e7;
          if(!config_setting_lookup_float(bunch_setting, "gamma0_std", &gamma0_std))
              gamma0_std = 0.0;
          
          // Spin parameters
          if(!config_setting_lookup_float(bunch_setting, "theta0_mean", &theta0_mean))
              theta0_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "theta0_std", &theta0_std))
              theta0_std = 0.0;
          
          if(!config_setting_lookup_float(bunch_setting, "phi0_mean", &phi0_mean))
              phi0_mean = 0.0;
          if(!config_setting_lookup_float(bunch_setting, "phi0_std", &phi0_std))
              phi0_std = 0.0;
          
          // Convert angles from degrees to radians for mean values
          theta0_mean *= M_PI/180.0;
          phi0_mean *= M_PI/180.0;
          theta0_std *= M_PI/180.0;  // Standard deviation also in radians
          phi0_std *= M_PI/180.0;
          
          printf("Bunch parameters:\n");
          printf("  zeta: mean=%f, std=%f\n", zeta0_mean, zeta0_std);
          printf("  x: mean=%e, std=%e\n", x_mean, x_std);
          printf("  y: mean=%e, std=%e\n", y_mean, y_std);
          printf("  px: mean=%e, std=%e\n", px_mean, px_std);
          printf("  py: mean=%e, std=%e\n", py_mean, py_std);
          printf("  gamma0: mean=%e, std=%e\n", gamma0_mean, gamma0_std);
          printf("  theta0: mean=%f, std=%f\n", theta0_mean, theta0_std);
          printf("  phi0: mean=%f, std=%f\n", phi0_mean, phi0_std);
          
          // Use reference gamma for time step calculation
          double ref_gamma0 = gamma0_mean;
          double r_omega_beta = sqrt(ref_gamma0/kappa_square);
          double dt = 2*M_PI*r_omega_beta/sampling_factor;
          fprintf(stdout,"dt = %f\n", dt);
          
          int total_time_steps = ((int) (t_duration / dt)) + 1;
          int storage_time_steps = ((int) (total_time_steps / storage_sampling_factor)) + 1;
          
          fprintf(stdout,"total_time_steps = %d\n", total_time_steps);
          fprintf(stdout,"storage_time_steps = %d (reduction factor: %.1f)\n", 
                  storage_time_steps, (double)total_time_steps/storage_time_steps);
          
          int column_length = 10; // t, zeta, x, y, pz, px, py, sz, sx, sy
          
          // Create time array for storage points
          double *time_array = (double*)malloc(storage_time_steps * sizeof(double));
          if (time_array == NULL)
          {
              fprintf(stderr, "Failed to allocate memory for time array\n");
              return EXIT_FAILURE;
          }
          
          for(int j = 0; j < storage_time_steps; j++)
          {
              time_array[j] = t0 + j * dt * storage_sampling_factor;
          }
          
          // Save time array to HDF5
          hsize_t time_dims[1] = {storage_time_steps};
          h5status = H5LTmake_dataset_double(ofile_h5id, "time", 1, time_dims, time_array);
          free(time_array);
          
          // Set batch size for memory management
          int batch_size = 50;
          int num_batches = (bunch_size + batch_size - 1) / batch_size; // Ceiling division
          
          printf("Processing %d particles in %d batches (batch size: %d)\n", 
                bunch_size, num_batches, batch_size);
          
          // Create the main dataset with the storage dimensions
          hsize_t full_dims[3] = {bunch_size, storage_time_steps, column_length};
          hid_t file_space = H5Screate_simple(3, full_dims, NULL);
          hid_t dset = H5Dcreate2(ofile_h5id, "bunch_data", H5T_NATIVE_DOUBLE, 
                                file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          
          // Save bunch parameters as attributes
          H5LTset_attribute_int(ofile_h5id, "bunch_data", "bunch_size", &bunch_size, 1);
          H5LTset_attribute_double(ofile_h5id, "bunch_data", "zeta0_mean", &zeta0_mean, 1);
          H5LTset_attribute_double(ofile_h5id, "bunch_data", "zeta0_std", &zeta0_std, 1);
          H5LTset_attribute_double(ofile_h5id, "bunch_data", "gamma0_mean", &gamma0_mean, 1);
          H5LTset_attribute_string(ofile_h5id, "bunch_data", "distribution", distribution_type);
          H5LTset_attribute_double(ofile_h5id, "bunch_data", "storage_sampling_factor", &storage_sampling_factor, 1);
          
          // Save variable names as attributes
          const char *var_names[] = {"t", "zeta", "x", "y", "pz", "px", "py", "sz", "sx", "sy"};
          for (int v = 0; v < 10; v++) 
          {
              char attr_name[50];
              sprintf(attr_name, "variable_name_%d", v);
              H5LTset_attribute_string(ofile_h5id, "bunch_data", attr_name, var_names[v]);
          }
          
          // Process particles in batches
          for(int batch = 0; batch < num_batches; batch++)
          {
              int batch_start = batch * batch_size;
              int batch_end = (batch + 1) * batch_size;
              if(batch_end > bunch_size) batch_end = bunch_size;
              int current_batch_size = batch_end - batch_start;
              
              printf("Processing batch %d/%d (particles %d-%d)\n", 
                    batch+1, num_batches, batch_start, batch_end-1);
              
              // Allocate memory only for current batch (storage resolution)
              double *batch_buffer = (double*)malloc(current_batch_size * storage_time_steps * column_length * sizeof(double));
              if (batch_buffer == NULL)
              {
                  fprintf(stderr, "Failed to allocate memory for batch %d\n", batch);
                  H5Dclose(dset);
                  H5Sclose(file_space);
                  return EXIT_FAILURE;
              }
              
              // Track particles in current batch using OpenMP
              #pragma omp parallel for schedule(dynamic)
              for(int local_i = 0; local_i < current_batch_size; local_i++)
              {
                  int global_i = batch_start + local_i;
                  unsigned int thread_seed = seed + global_i;
                  
                  // Generate random initial conditions based on distribution type
                  double zeta0, x0, y0, px0, py0, gamma0, theta0, phi0;
                  
                  if(strcmp(distribution_type, "gaussian") == 0)
                  {
                      zeta0 = gaussian_rand_threadsafe(zeta0_mean, zeta0_std, &thread_seed);
                      x0 = gaussian_rand_threadsafe(x_mean, x_std, &thread_seed);
                      y0 = gaussian_rand_threadsafe(y_mean, y_std, &thread_seed);
                      px0 = gaussian_rand_threadsafe(px_mean, px_std, &thread_seed);
                      py0 = gaussian_rand_threadsafe(py_mean, py_std, &thread_seed);
                      gamma0 = gaussian_rand_threadsafe(gamma0_mean, gamma0_std, &thread_seed);
                      theta0 = gaussian_rand_threadsafe(theta0_mean, theta0_std, &thread_seed);
                      phi0 = gaussian_rand_threadsafe(phi0_mean, phi0_std, &thread_seed);
                  }
                  else // uniform distribution
                  {
                      zeta0 = uniform_rand_threadsafe(zeta0_mean - zeta0_std, zeta0_mean + zeta0_std, &thread_seed);
                      x0 = uniform_rand_threadsafe(x_mean - x_std, x_mean + x_std, &thread_seed);
                      y0 = uniform_rand_threadsafe(y_mean - y_std, y_mean + y_std, &thread_seed);
                      px0 = uniform_rand_threadsafe(px_mean - px_std, px_mean + px_std, &thread_seed);
                      py0 = uniform_rand_threadsafe(py_mean - py_std, py_mean + py_std, &thread_seed);
                      gamma0 = uniform_rand_threadsafe(gamma0_mean - gamma0_std, gamma0_mean + gamma0_std, &thread_seed);
                      theta0 = uniform_rand_threadsafe(theta0_mean - theta0_std, theta0_mean + theta0_std, &thread_seed);
                      phi0 = uniform_rand_threadsafe(phi0_mean - phi0_std, phi0_mean + phi0_std, &thread_seed);
                  }
                  
                  // Ensure gamma0 is positive
                  if(gamma0 < 1.0) gamma0 = 1.0;
                  
                  // Add pz0 protection
                  double energy_term = gamma0*gamma0 - 1.0 - px0*px0 - py0*py0;
                  if(energy_term < 0) {
                      // If energy term is negative, adjust transverse momenta
                      double max_transverse_momentum = sqrt(gamma0*gamma0 - 1.0) * 0.9; // Reserve 10% for longitudinal
                      double current_transverse = sqrt(px0*px0 + py0*py0);
                      if(current_transverse > max_transverse_momentum) {
                          double scale = max_transverse_momentum / current_transverse;
                          px0 *= scale;
                          py0 *= scale;
                      }
                  }
                  
                  // Calculate initial pz from energy-momentum relation
                  double pz0 = sqrt(gamma0*gamma0 - 1.0 - px0*px0 - py0*py0);
                  
                  // Set initial conditions in batch buffer for this particle
                  int particle_offset = local_i * storage_time_steps * column_length;
                  
                  // First storage point (j=0)
                  int storage_index = 0;
                  int storage_offset = particle_offset + storage_index * column_length;
                  batch_buffer[storage_offset + 0] = t0;           // t
                  batch_buffer[storage_offset + 1] = zeta0;        // zeta
                  batch_buffer[storage_offset + 2] = x0;           // x
                  batch_buffer[storage_offset + 3] = y0;           // y
                  batch_buffer[storage_offset + 4] = pz0;          // pz
                  batch_buffer[storage_offset + 5] = px0;          // px
                  batch_buffer[storage_offset + 6] = py0;          // py
                  batch_buffer[storage_offset + 7] = cos(theta0);  // sz
                  batch_buffer[storage_offset + 8] = sin(theta0) * cos(phi0);  // sx
                  batch_buffer[storage_offset + 9] = sin(theta0) * sin(phi0);  // sy
                  
                  // Create temporary arrays for full-resolution tracking
                  double *current_state = (double*)malloc(column_length * sizeof(double));
                  double *next_state = (double*)malloc(column_length * sizeof(double));
                  
                  // Initialize current state
                  for(int k = 0; k < column_length; k++) {
                      current_state[k] = batch_buffer[storage_offset + k];
                  }
                  
                  // Track the particle with full resolution but store only at intervals
                  double current_time = t0;
                  int storage_counter = 1; // Next storage index
                  
                  for(int step = 1; step < total_time_steps && storage_counter < storage_time_steps; step++)
                  {
                      double next_time = current_time + dt;
                      
                      // Perform RK4 step
                      int opt = rk4vec(current_time, column_length-1, 
                                      current_state + 1, 
                                      next_state + 1, dt, dy_over_dt);
                      
                      if(opt == -1)
                      {
                          // Error handling - copy previous state
                          for(int k = 1; k < column_length; k++) {
                              next_state[k] = current_state[k];
                          }
                      }
                      
                      next_state[0] = next_time;
                      
                      // Check if we should store this step
                      if(fmod(step, storage_sampling_factor) < 1.0) 
                      {
                          int new_storage_offset = particle_offset + storage_counter * column_length;
                          for(int k = 0; k < column_length; k++) {
                              batch_buffer[new_storage_offset + k] = next_state[k];
                          }
                          storage_counter++;
                      }
                      
                      // Update current state for next iteration
                      for(int k = 0; k < column_length; k++) {
                          current_state[k] = next_state[k];
                      }
                      current_time = next_time;
                  }
                  
                  // Ensure we store the final state if we haven't already
                  if(storage_counter < storage_time_steps) {
                      int final_storage_offset = particle_offset + storage_counter * column_length;
                      for(int k = 0; k < column_length; k++) {
                          batch_buffer[final_storage_offset + k] = current_state[k];
                      }
                  }
                  
                  free(current_state);
                  free(next_state);
              }
              
              // Write current batch to HDF5 file
              hsize_t start[3] = {batch_start, 0, 0};
              hsize_t count[3] = {current_batch_size, storage_time_steps, column_length};
              
              // Select hyperslab in the file
              H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
              
              // Create memory space for the batch
              hsize_t mem_dims[3] = {current_batch_size, storage_time_steps, column_length};
              hid_t mem_space = H5Screate_simple(3, mem_dims, NULL);
              
              // Write batch data to file
              h5status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, 
                                H5P_DEFAULT, batch_buffer);
              
              if(h5status < 0) {
                  fprintf(stderr, "Error writing batch %d to HDF5 file\n", batch);
              }
              
              // Clean up memory space and batch buffer
              H5Sclose(mem_space);
              free(batch_buffer);
              
              printf("Batch %d/%d completed and saved\n", batch+1, num_batches);
              
              // Optional: Flush HDF5 file to ensure data is written to disk
              H5Fflush(ofile_h5id, H5F_SCOPE_GLOBAL);
          }
          
          // Close HDF5 resources
          H5Dclose(dset);
          H5Sclose(file_space);
          
          printf("All %d batches processed successfully. Data saved to HDF5 file.\n", num_batches);
          printf("Storage reduction: from %d to %d time steps (factor: %.1f)\n", 
                total_time_steps, storage_time_steps, (double)total_time_steps/storage_time_steps);
      }
      break;

  }

  config_destroy(&cfg);
  h5status = H5Fclose (ofile_h5id);
  printf("h5status = %d\n", h5status);
  return(EXIT_SUCCESS);
}

/* eof */