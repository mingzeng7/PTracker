PTracker is a particle tracker in a blowout regime plasma wakefield accelerator. The units are normalized: length is normalized to k_p^-1, time is normalized to omega_p^-1, force is normalized to m_e*c^2*k_p, momentum is normalized to m_e*c.

To compile, one has to install hdf5 (recommend version 1.10) and libconfig, and set their lib and header file path in Makefile, and then use "make" command. Then a excutable "PTracker" will be created. It uses "input.cfg" as the default input deck. 

The code is solving the equations f_x_ext = -x/2, f_z_ext = -zeta/2, where zeta = z - beta_w*t, beta_w is the phase speed of the wake normalized to c. One can specify wheter to turn on radiation reaction by setting "if_RR" in the input deck.

Currently only single particle mode is implemented.

One can use the Python script "plotter_PTracker.py" for plotting the trajectory.
