# Display spin precession characteristics

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import find_peaks
import pylibconfig2 as libconfig
from datetime import datetime
from scipy.integrate import cumtrapz

class PT_plt:
    def __init__(self, h_file, particle_name='0', config=None):
        self.h_file = h_file
        self.particle_name = particle_name
        self.axis_dict = {'t':0, 'zeta':1, 'x':2, 'y':3, 'pz':4, 'px':5, 'py':6, 'sz':7, 'sx':8, 'sy':9}
        self.data = self.h_file[self.particle_name]
        
        # Store configuration parameters
        if config:
            self.kappa_square = config.get('kappa_square', 0.5)
            self.lambd = config.get('lambd', 0.0)
            self.f_z0 = config.get('f_z0', 0.0)
            self.elec_anomaly = config.get('elec_anomaly', 1.1596521884e-3)
            self.sampling_factor = config.get('sampling_factor', 1.0)
        else:
            self.kappa_square = 0.5
            self.lambd = 0.0
            self.f_z0 = 0.0
            self.elec_anomaly = 1.1596521884e-3
            self.sampling_factor = 1.0
        
        # Initialize other attributes
        self.gama = None
        self.indices = None
        self.omega_beta = None
        
    def get(self, axis):
        return self.data[:, self.axis_dict[axis]]
        
    def get_half_step(self, axis):
        return (self.data[1:, self.axis_dict[axis]]+self.data[:-1, self.axis_dict[axis]])/2
        
    def get_d_half_step(self, axis):
        return self.data[1:, self.axis_dict[axis]]-self.data[:-1, self.axis_dict[axis]]
        
    def save_gamma(self):
        self.gama = np.sqrt(np.square(self.get('pz')) + np.square(self.get('px')) + np.square(self.get('py')) + 1.)
        
    def save_get_gamma(self):
        if self.gama is None:
            self.save_gamma()
        return self.gama
        
    def get_gamma_half_step(self):
        gama = self.save_get_gamma()
        return (gama[1:]+gama[:-1])/2
        
    def get_d_gamma_half_step(self):
        gama = self.save_get_gamma()
        return gama[1:]-gama[:-1]
        
    def find_peak_inds(self, axis):
        return find_peaks(np.abs(self.get(axis)))[0]
    
    def find_valley_inds(self, axis):
        return find_peaks(-self.get(axis))[0]
    
    def find_pediment_inds(self, axis):
        return find_peaks(self.get(axis))[0]
    
    def get_axis_average(self, axis):
        t = self.get('t')
        sx_pediment_inds = self.find_pediment_inds(axis)
        sx_valley_inds = self.find_valley_inds(axis)
        sx_valley = self.get(axis)[sx_valley_inds]
        sx_pediment = self.get(axis)[sx_pediment_inds]
        t_sx_pediment = t[sx_pediment_inds]
        sx_valley_interp = np.interp(t_sx_pediment, t[sx_valley_inds], sx_valley)
        return (sx_pediment+sx_valley_interp)/2, t_sx_pediment
    
    def get_spin_x_average(self):
        return self.get_axis_average('sx')
    
    def get_spin_y_average(self):
        return self.get_axis_average('sy')
    
    def get_spin_z_average(self):
        return self.get_axis_average('sz')
    
    def get_x1_times_px1(self):
        t = self.get('t')
        x_peak_inds = self.find_peak_inds('x')
        x_peaks = np.abs(self.get('x')[x_peak_inds])
        px_peak_inds = self.find_peak_inds('px')
        px_peaks = np.abs(self.get('px')[px_peak_inds])
        t_x_peaks = t[x_peak_inds]
        px_peaks_interp = np.interp(t_x_peaks, t[px_peak_inds], px_peaks)
        return x_peaks*px_peaks_interp, t_x_peaks
        
    def get_y1_times_py1(self):
        t = self.get('t')
        y_peak_inds = self.find_peak_inds('y')
        y_peaks = np.abs(self.get('y')[y_peak_inds])
        py_peak_inds = self.find_peak_inds('py')
        py_peaks = np.abs(self.get('py')[py_peak_inds])
        t_y_peaks = t[y_peak_inds]
        py_peaks_interp = np.interp(t_y_peaks, t[py_peak_inds], py_peaks)
        return y_peaks*py_peaks_interp, t_y_peaks
        
    def get_Sx(self, if_multi_pi=True):
        Sx, t = self.get_x1_times_px1()
        if if_multi_pi:
            Sx *= np.pi
        return Sx, t
        
    def get_Sy(self, if_multi_pi=True):
        Sy, t = self.get_y1_times_py1()
        if if_multi_pi:
            Sy *= np.pi
        return Sy, t
    
    def get_sx_peaks_T(self):
        t = self.get('t')
        sx_peak_inds = self.find_peak_inds('sx')
        sx_peaks = np.abs(self.get('sx')[sx_peak_inds])
        t_sx_peaks = t[sx_peak_inds]
        sx_peak_inds1 = find_peaks(sx_peaks)[0]
        sx_peaks1 = sx_peaks[sx_peak_inds1]
        t_sx_peaks1 = t_sx_peaks[sx_peak_inds1]
        return sx_peaks1, t_sx_peaks1
    
    def get_E_z(self):
        zeta = self.get('zeta')
        return -self.f_z0 + self.lambd * zeta
   
    def get_omega_beta(self):
        gama = self.save_get_gamma()
        return np.sqrt(self.kappa_square/gama)
    
    def get_U(self):
        gama = self.save_get_gamma()
        x, beta_x, t = self.get('x'), self.get('px')/gama, self.get('t')
        omega_beta = np.sqrt(self.kappa_square/gama)
        phi = cumtrapz(omega_beta, t, initial=0)
        return (x-(beta_x/omega_beta)*1j)*np.exp(-(phi)*1j)
    
    def get_V(self):
        gama = self.save_get_gamma()
        y, beta_y, t = self.get('y'), self.get('py')/gama, self.get('t')
        omega_beta = np.sqrt(self.kappa_square/gama)
        phi = cumtrapz(omega_beta, t, initial=0)
        return (y-(beta_y/omega_beta)*1j)*np.exp(-(phi)*1j)
    
    def get_lz(self):
        x = self.get('x')
        y = self.get('y')
        px = self.get('px')
        py = self.get('py')   
        return x*py - y*px
    
    def set_indices(self, indices):
        self.indices = indices
    
    def average_around_k(self, array, k):
        """Calculate average of array elements around index k within window_size radius"""
        if self.omega_beta is None:
            self.omega_beta = self.get_omega_beta()
        
        n = 1  # n refers to half window_size in number of periods
        window_size = int(self.sampling_factor * n * self.omega_beta[0] / self.omega_beta[k])
        start_index = k - window_size
        end_index = k + window_size + 1
        
        while start_index < 0 or end_index > len(array):
            n -= 1
            window_size = int(self.sampling_factor * n * self.omega_beta[0] / self.omega_beta[k])
            start_index = k - window_size
            end_index = k + window_size + 1

        sub_array = array[start_index:end_index]
        average = np.mean(sub_array)
        return average
    
    def average_function(self, array, indices=None):
        if indices is None:
            if self.indices is None:
                raise ValueError("Indices not set. Use set_indices() or provide indices parameter.")
            indices = self.indices
        
        sub_array = []
        for i in indices[0:]:
            sub_array += [self.average_around_k(array, i)]
        return np.array(sub_array)
    
    def get_Omega_T(self, indices=None):
        gama = self.save_get_gamma()
        pz = self.get('pz')
        beta_z = pz/gama
        E_z = self.get_E_z()
        lz = self.get_lz()
        
        gama_average = self.average_function(gama, indices)
        beta_z0 = self.average_function(beta_z, indices)
        E_z0 = self.average_function(E_z, indices)
        lz_average = self.average_function(lz, indices)
        
        Omega_T = -(self.kappa_square*(self.elec_anomaly+(1-self.lambd)/gama_average)/gama_average - 
                   (gama_average*beta_z0**2*self.kappa_square+E_z0**2)*(self.elec_anomaly+1/gama_average)**2/2/gama_average)*lz_average
        return Omega_T
    
    def get_Omega_T_direct_average_method(self):
        t = self.get('t')
        gama = self.save_get_gamma()
        pz = self.get('pz')
        beta_z = pz/gama
        E_z = self.get_E_z()
        lz = self.get_lz()
        
        beta_z0, t_beta_z0 = find_array_average(t, beta_z)
        E_z0, t_Ez0 = find_array_average(t, E_z)
        lz0, t_lz0 = find_array_average(t, lz)
        beta_z0_interp = np.interp(t, t_beta_z0, beta_z0)
        E_z0_interp = np.interp(t, t_Ez0, E_z0)
        lz0_interp = np.interp(t, t_lz0, lz0)
        
        if self.indices is None:
            raise ValueError("Indices not set. Use set_indices() first.")
        
        gama_average = gama[self.indices]
        beta_z0 = beta_z0_interp[self.indices]
        E_z0 = E_z0_interp[self.indices]
        lz_average = lz0_interp[self.indices]
        
        Omega_T = -(self.kappa_square*(self.elec_anomaly+(1-self.lambd)/gama_average)/gama_average + 
                   (gama_average*beta_z0**2*self.kappa_square+E_z0**2)*(self.elec_anomaly+1/gama_average)**2/2/gama_average)*lz_average
        return Omega_T
    
    def get_theory_s1x(self, s0z, indices=None):
        t = self.get('t')
        gama = self.save_get_gamma()
        pz = self.get('pz')
        px = self.get('px')
        beta_z = pz/gama
        beta_x = px/gama
        E_z = self.get_E_z()
        
        if indices is None:
            if self.indices is None:
                raise ValueError("Indices not set. Use set_indices() or provide indices parameter.")
            indices = self.indices
        
        gama_average = self.average_function(gama, indices)
        beta_z0 = self.average_function(beta_z, indices)
        E_z0 = self.average_function(E_z, indices)
        
        T = t[indices]
        gama_average_interp = np.interp(t, T, gama_average)
        beta_z0_interp = np.interp(t, T, beta_z0)
        E_z0_interp = np.interp(t, T, E_z0)
        
        x = self.get('x')
        return (self.elec_anomaly+1/gama_average_interp)*s0z*(gama_average_interp*beta_z0_interp*beta_x+E_z0_interp*x)
    
    def get_theory_s1y(self, s0z, indices=None):
        t = self.get('t')
        gama = self.save_get_gamma()
        pz = self.get('pz')
        py = self.get('py')
        beta_z = pz/gama
        beta_y = py/gama
        E_z = self.get_E_z()
        
        if indices is None:
            if self.indices is None:
                raise ValueError("Indices not set. Use set_indices() or provide indices parameter.")
            indices = self.indices
        
        gama_average = self.average_function(gama, indices)
        beta_z0 = self.average_function(beta_z, indices)
        E_z0 = self.average_function(E_z, indices)
        
        T = t[indices]
        gama_average_interp = np.interp(t, T, gama_average)
        beta_z0_interp = np.interp(t, T, beta_z0)
        E_z0_interp = np.interp(t, T, E_z0)
        
        y = self.get('y')
        return (self.elec_anomaly+1/gama_average_interp)*s0z*(gama_average_interp*beta_z0_interp*beta_y+E_z0_interp*y)
    
    def get_theory_s1z(self, s0x, s0y, indices=None):
        t = self.get('t')
        gama = self.save_get_gamma()
        pz = self.get('pz')
        px = self.get('px')
        py = self.get('py')
        beta_z = pz/gama
        beta_x = px/gama
        beta_y = py/gama
        E_z = self.get_E_z()
        
        if indices is None:
            if self.indices is None:
                raise ValueError("Indices not set. Use set_indices() or provide indices parameter.")
            indices = self.indices
        
        gama_average = self.average_function(gama, indices)
        beta_z0 = self.average_function(beta_z, indices)
        E_z0 = self.average_function(E_z, indices)
        
        T = t[indices]
        gama_average_interp = np.interp(t, T, gama_average)
        beta_z0_interp = np.interp(t, T, beta_z0)
        E_z0_interp = np.interp(t, T, E_z0)
        
        y = self.get('y')
        x = self.get('x')
        return -(self.elec_anomaly+1/gama_average_interp)*(s0y*(gama_average_interp*beta_z0_interp*beta_y+E_z0_interp*y)+
                                                          s0x*(gama_average_interp*beta_z0_interp*beta_x+E_z0_interp*x))
    
    def get_sin_Delta_Phi(self, indices=None):
        U = self.get_U()
        V = self.get_V()
        lz = self.get_lz()
        gama = self.save_get_gamma()
        
        if indices is None:
            if self.indices is None:
                raise ValueError("Indices not set. Use set_indices() or provide indices parameter.")
            indices = self.indices
        
        lz_average = self.average_function(lz, indices)
        U_average = self.average_function(U, indices)
        V_average = self.average_function(V, indices)
        gama_average = self.average_function(gama, indices)
        
        U_abs = np.abs(U_average)
        V_abs = np.abs(V_average)
        sin_Delta_Phi = -lz_average/np.sqrt(self.kappa_square*gama_average)/U_abs/V_abs
        
        return sin_Delta_Phi


def caculate_s0(s0, Phi):
    sx = s0[0]*np.cos(Phi) + s0[1]*np.sin(Phi)
    sy = -s0[0]*np.sin(Phi) + s0[1]*np.cos(Phi)
    sz = s0[2]* np.cos(Phi-Phi)
    return sx, sy, sz


def find_array_average(t, array):
    sx_pediment_inds = find_peaks(array)[0]
    sx_valley_inds = find_peaks(-array)[0]
    sx_valley = array[sx_valley_inds]
    sx_pediment = array[sx_pediment_inds]
    t_sx_pediment = t[sx_pediment_inds]
    sx_valley_interp = np.interp(t_sx_pediment, t[sx_valley_inds], sx_valley)
    return (sx_pediment+sx_valley_interp)/2, t_sx_pediment


def particle_trace_long_term_equation(t, y, r_e, kappa_square, gama_w, f_z0, lambd):
    y1, y2, y3, y4 = y
    dy1_dt = 1/2/gama_w**2 - 1/2/y4**2 - 1/4*np.sqrt(kappa_square/y4**3)*y2
    dy2_dt = -1/4*r_e*np.sqrt(kappa_square**3*y4)*(y2**2+4/3*y3**2)
    dy3_dt = -1/3*r_e*np.sqrt(kappa_square**3*y4)*y2*y3
    dy4_dt = -(-f_z0+lambd*y1)*(1-1/2/y4**2-1/4*np.sqrt(kappa_square/y4**3)*y2) - 1/3*r_e*np.sqrt(kappa_square**3*y4**3)*y2
    return np.array([dy1_dt, dy2_dt, dy3_dt, dy4_dt])


def rk4_system(f, y0, t, *args):
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0

    for i in range(n - 1):
        h = t[i + 1] - t[i]
        k1 = f(t[i], y[i], *args)
        k2 = f(t[i] + h / 2, y[i] + h / 2 * k1, *args)
        k3 = f(t[i] + h / 2, y[i] + h / 2 * k2, *args)
        k4 = f(t[i] + h, y[i] + h * k3, *args)
        y[i + 1] = y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return y


def caculate_Omega_T(y_rk, kappa_square, elec_anomaly, lambd, f_z0):
    zeta_average = y_rk[:, 0]
    S_average = y_rk[:, 1]
    lz_average = y_rk[:, 2]
    gama_average = y_rk[:, 3]
    beta_z0 = 1 - 1/2/gama_average**2 - 1/4*np.sqrt(kappa_square/gama_average**6)*S_average
    E_z0 = -f_z0 + lambd*zeta_average
    Omega_T = -(kappa_square*(elec_anomaly+(1-lambd)/gama_average)/gama_average - 
               (gama_average*beta_z0**2*kappa_square+E_z0**2)*(elec_anomaly+1/gama_average)**2/2/gama_average)*lz_average
    return Omega_T


def IF_RR(if_RR):
    return 1.0 if if_RR else 0.0


def improved_fft_analysis(signal, dt, window='hann', detrend=True):
    n = len(signal)
    
    signal_amplitude = np.max(np.abs(signal))
    signal_mean = np.mean(signal)
    
    if detrend:
        signal_detrended = signal - np.mean(signal)
    else:
        signal_detrended = signal
    
    if window == 'hann':
        window_func = np.hanning(n)
    elif window == 'hamming':
        window_func = np.hamming(n)
    elif window == 'blackman':
        window_func = np.blackman(n)
    else:
        window_func = np.ones(n)
    
    window_energy = np.sum(window_func**2)
    window_coherent_gain = np.sum(window_func) / n
    
    signal_windowed = signal_detrended * window_func
    
    fft_result = np.fft.fft(signal_windowed)
    fft_amp_standard = 2.0 / (n * window_coherent_gain) * np.abs(fft_result[:n//2])
    freqs = np.fft.fftfreq(n, dt)[:n//2]
    
    return freqs, fft_amp_standard


def find_dominant_frequency_range(freqs, fft_amp, threshold_ratio=0.1):
    peak_idx = np.argmax(fft_amp)
    peak_freq = freqs[peak_idx]
    peak_amp = fft_amp[peak_idx]
    
    threshold = peak_amp * threshold_ratio
    significant_indices = np.where(fft_amp >= threshold)[0]
    
    if len(significant_indices) == 0:
        return (freqs[0], freqs[-1]), peak_freq
    
    min_freq_idx = np.min(significant_indices)
    max_freq_idx = np.max(significant_indices)
    
    min_freq = freqs[min_freq_idx]
    max_freq = freqs[max_freq_idx]
    
    return (min_freq, max_freq), peak_freq


def get_optimal_freq_range(peak_freq_theory, peak_freq_nume, scale_factor=5.0):
    max_peak_freq = max(peak_freq_theory, peak_freq_nume)
    freq_min = 0
    freq_max = max_peak_freq * scale_factor
    
    if freq_max < 0.001:
        freq_max = 0.01
    
    return (freq_min, freq_max)


if __name__ == "__main__":
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    config = libconfig.Config()
    config.read_file('input.cfg')
    
    output_file = config.lookup('output_file')
    if_RR = config.lookup('if_RR')
    if_RR1 = config.lookup('if_RR1')
    if_SG = config.lookup('if_SG')
    if_SG_grad = config.lookup('if_SG_grad')
    
    r_e = 2.8179403262e-15 * 1e5 * IF_RR(if_RR)
    elec_anomaly = 1.1596521884e-3
    sampling_factor = config.lookup('sampling_factor')
    t1 = config.lookup('t1')
    
    particle = config.lookup('particle')
    kappa_square = 0.5
    
    single_particle = particle._lookup(['single_particle'])[0]
    U = single_particle._lookup(['x_1'])
    V = single_particle._lookup(['y_1'])
    gama0 = single_particle._lookup(['gamma0'])
    zeta0 = single_particle._lookup(['zeta0'])
    phase0_x = single_particle._lookup(['phase0_x'])
    phase0_y = single_particle._lookup(['phase0_y'])
    theta0 = single_particle._lookup(['theta0'])
    phi0 = single_particle._lookup(['phi0'])
    
    wake = config.lookup('wake')
    kp_SI = wake._lookup(['kp_SI'])
    gama_w = wake._lookup(['gamma_w'])
    lambd = wake._lookup(['lambda'])
    f_z0 = wake._lookup(['f_z0'])
    
    # Create config dictionary for PT_plt
    pt_config = {
        'kappa_square': kappa_square,
        'lambd': lambd,
        'f_z0': f_z0,
        'elec_anomaly': elec_anomaly,
        'sampling_factor': sampling_factor
    }
    
    # Save config backup
    try:
        with open('input.cfg', 'r') as cfg_file:
            cfg_content = cfg_file.read()
        
        cfg_backup_filename = f"{output_file.replace('.h5', '')}_config_backup_{timestamp}.txt"
        
        with open(cfg_backup_filename, 'w') as backup_file:
            backup_file.write(f"Configuration Backup - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            backup_file.write("=" * 50 + "\n\n")
            backup_file.write(cfg_content)
        
        print(f"Configuration backup saved as: {cfg_backup_filename}")
    except Exception as e:
        print(f"Warning: Could not backup configuration file. Error: {e}")
    
    omega_beta_0 = np.sqrt(kappa_square/gama0)
    Omega_T0 = (kappa_square*(elec_anomaly+(1-lambd)/gama0)/gama0 + 
               (gama0*kappa_square+(-f_z0 + lambd * zeta0)**2)*(elec_anomaly+1/gama0)**2/2/gama0) * \
               np.sqrt(kappa_square*gama0*U*V)*np.sin((phase0_y-phase0_x)/180*np.pi)
    
    # Initial conditions for particle long-term equations
    y0 = np.array([zeta0, np.sqrt(kappa_square*gama0)*(U**2+V**2), 
                   -np.sqrt(kappa_square*gama0)*U*V*np.sin((phase0_y-phase0_x)*np.pi/180), gama0])
    t_trace = np.linspace(0, t1, int(t1/100))
    
    # Solve long-term equations with parameters passed as arguments
    y_rk4 = rk4_system(particle_trace_long_term_equation, y0, t_trace, r_e, kappa_square, gama_w, f_z0, lambd)
    
    # Calculate Omega_T from long-term equations
    Omega_T_trace_method = caculate_Omega_T(y_rk4, kappa_square, elec_anomaly, lambd, f_z0)
    
    with h5py.File(output_file, 'r') as h_file:
        # Create PT_plt instance with config
        plotter = PT_plt(h_file, config=pt_config)
        gama = plotter.save_get_gamma()
        
        sx, t = plotter.get('sx'), plotter.get('t')
        sz = plotter.get('sz')
        sy = plotter.get('sy')
        
        Dt = np.abs(2*np.pi/Omega_T0/sampling_factor)
        dt = t[1] - t[0]
        n = int(Dt/dt)
        indices = np.linspace(0, len(t)-1, n, dtype=int)
        
        # Set indices for the plotter
        plotter.set_indices(indices)
        
        # Calculate Omega_T using the plotter
        Omega_T = plotter.get_Omega_T()
        
        print(f"Omega_T0: {Omega_T0}, Omega_T[0]: {Omega_T[0]}")
        
        T = t[indices]
        print(f"T[0]: {T[0]}")
        
        Phi = cumtrapz(Omega_T, T, initial=0)
        Phi_trace_method = cumtrapz(Omega_T_trace_method, t_trace, initial=0)
        
        # Calculate spin averages
        sx_average, t_x_average = plotter.get_spin_x_average()
        sy_average, t_y_average = plotter.get_spin_y_average()
        sz_average, t_z_average = plotter.get_spin_z_average()
        
        sx_average_interp = np.interp(t, t_x_average, sx_average)
        sy_average_interp = np.interp(t, t_y_average, sy_average)
        sz_average_interp = np.interp(t, t_z_average, sz_average)
        
        # Calculate theory vs numerical spin components
        s1x_theory = plotter.get_theory_s1x(sz_average_interp)
        s1x_nume = sx - sx_average_interp
        s1x_nume_peak_inds = find_peaks(s1x_nume)[0]
        s1x_nume_peaks = s1x_nume[s1x_nume_peak_inds]
        t_x_nume_peaks = t[s1x_nume_peak_inds]
        s1y_theory = plotter.get_theory_s1y(sz_average_interp)
        s1y_nume = sy - sy_average_interp
        s1z_theory = plotter.get_theory_s1z(sx_average_interp, sy_average_interp)
        s1z_nume = sz - sz_average_interp
        s1z_nume_peak_inds = find_peaks(s1z_nume)[0]
        s1z_nume_peaks = s1z_nume[s1z_nume_peak_inds]
        t_z_nume_peaks = t[s1z_nume_peak_inds]

        # Calculate spin from long-term theory
        sx0_tm, sy0_tm, sz0_tm = caculate_s0([sx[0]-s1x_theory[0], sy[0]-s1y_theory[0], sz_average_interp[0]], Phi_trace_method)
        sx0, sy0, sz0 = caculate_s0([sx_average_interp[0], sy_average_interp[0], sz_average_interp[0]], Phi)
        
        sx0_tm_interp = np.interp(t, t_trace, sx0_tm)
        sy0_tm_interp = np.interp(t, t_trace, sy0_tm)
        sz0_tm_interp = np.interp(t, t_trace, sz0_tm)
        
        s1x_theory_tm = plotter.get_theory_s1x(sz0_tm_interp)
        s1x_theory_peaks_inds = find_peaks(s1x_theory_tm)[0]
        s1x_theory_peaks = s1x_theory_tm[s1x_theory_peaks_inds]
        t_x_theory_peaks = t[s1x_theory_peaks_inds]
        s1y_theory_tm = plotter.get_theory_s1y(sz0_tm_interp)
        s1z_theory_tm = plotter.get_theory_s1z(sx0_tm_interp, sy0_tm_interp)
        s1z_theory_peaks_inds = find_peaks(s1z_theory_tm)[0]
        s1z_theory_peaks = s1z_theory_tm[s1z_theory_peaks_inds]
        t_z_theory_peaks = t[s1z_theory_peaks_inds]        
        print(f'Initial sx0 value: {sx_average_interp[0]:.6f}, Initial sy0 value: {sy_average_interp[0]:.6f}')
        
        sin_Delta_Phi = plotter.get_sin_Delta_Phi()
        
        
        plt.rcParams.update({
            'font.size': 28,  
            'axes.titlesize': 32,  
            'axes.labelsize': 32,  
            'xtick.labelsize': 28,  
            'ytick.labelsize': 28,  
            'legend.fontsize': 28,  
            'figure.titlesize': 34,  
            'lines.linewidth': 3.0,  
            'lines.markersize': 12,  
            'axes.linewidth': 2.5,  
            'grid.linewidth': 1.2,  
            'mathtext.default': 'regular', 
            'mathtext.fontset': 'stix',  
            'font.family': 'serif',
            'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
            'savefig.dpi': 600, 
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.02, 
            'figure.constrained_layout.use': False,  
        })

        colors = {
            'blue': '#1f77b4',
            'orange': '#ff7f0e',
            'green': '#2ca02c',
            'red': '#d62728',
            'purple': '#9467bd',
            'brown': '#8c564b',
            'pink': '#e377c2',
            'gray': '#7f7f7f',
            'yellow': '#bcbd22',
            'cyan': '#17becf'
        }

        # 1. comparison for s1z
        fig_s1z_time, ax_time_s1z = plt.subplots(figsize=(10, 7))  

        fig_width, fig_height = fig_s1z_time.get_size_inches()



        font_pt = plt.rcParams['font.size']
        label_pt = plt.rcParams['axes.labelsize']



        line2_s1z, = ax_time_s1z.plot(t, s1z_nume, color=colors['blue'], linewidth=1.2, alpha=0.5, label='Numerical')
        ax_time_s1z.set_xlabel('$t$', fontsize=50, labelpad=8)  
        ax_time_s1z.set_ylabel('$s_{1_z}$', fontsize=50, labelpad=8)  


        ax_time_s1z.text(0.02, 0.98, '(c)', transform=ax_time_s1z.transAxes, fontsize=42, fontweight='bold', 
                        verticalalignment='top', horizontalalignment='left')

        ax_time_s1z.grid(True, alpha=0.3, linestyle='--')
        ax_time_s1z.tick_params(axis='both', which='major', width=2.0, length=10, labelsize=40, pad=6)  
        ax_time_s1z.set_position([0.08, 0.10, 0.88, 0.85])  
        ax_time_s1z.legend([ line2_s1z], [ 'Numerical'],
                        loc='upper right', fontsize=40, frameon=True, framealpha=0.5,
                        edgecolor='black', fancybox=True,
                        bbox_to_anchor=(0.98, 0.92), borderaxespad=0.3)

        plt.subplots_adjust(left=0.06, right=0.98, top=0.98, bottom=0.08, hspace=0.01, wspace=0.01)

        import os
        output_dir, output_filename = os.path.split(output_file)
        plot_dir = output_dir + "single_particle_plot"
        os.makedirs(plot_dir, exist_ok=True)
        base_filename = os.path.splitext(output_filename)[0]
        s1z_time_base = f"{base_filename}_s1z_time_domain"
        s1z_time_path_png = os.path.join(plot_dir, f"{s1z_time_base}.png")
        s1z_time_path_pdf = os.path.join(plot_dir, f"{s1z_time_base}.pdf")
        plt.savefig(s1z_time_path_png, dpi=600, bbox_inches="tight", pad_inches=0.02)  
        plt.savefig(s1z_time_path_pdf, dpi=600, bbox_inches="tight", pad_inches=0.02)  
        plt.close(fig_s1z_time)

        # 2. comparison for s1x
        fig_s1x_time, ax_time_s1x = plt.subplots(figsize=(10, 7))

        line1_s1x, = ax_time_s1x.plot(t, s1x_theory_tm, color=colors['red'], linewidth=1.2, label='Theory')
        line2_s1x, = ax_time_s1x.plot(t, s1x_nume, color=colors['blue'], linewidth=1.2, alpha=0.5, label='Numerical')
        ax_time_s1x.set_xlabel('$t$', fontsize=50, labelpad=8)
        ax_time_s1x.set_ylabel('$s_{1_x}$', fontsize=50, labelpad=8)

        ax_time_s1x.text(0.02, 0.98, '(d)', transform=ax_time_s1x.transAxes, fontsize=40, fontweight='bold', 
                        verticalalignment='top', horizontalalignment='left')

        ax_time_s1x.grid(True, alpha=0.3, linestyle='--')
        ax_time_s1x.tick_params(axis='both', which='major', width=2.0, length=10, labelsize=40, pad=6)


        ax_time_s1x.set_position([0.08, 0.10, 0.88, 0.85])


        ax_time_s1x.legend([line1_s1x, line2_s1x], ['Theory', 'Numerical'],
                        loc='upper right', fontsize=40, frameon=True, framealpha=0.5,
                        edgecolor='black', fancybox=True,
                        bbox_to_anchor=(0.98, 0.92), borderaxespad=0.3)

        plt.subplots_adjust(left=0.06, right=0.98, top=0.98, bottom=0.08, hspace=0.01, wspace=0.01)

        s1x_time_base = f"{base_filename}_s1x_time_domain"
        s1x_time_path_png = os.path.join(plot_dir, f"{s1x_time_base}.png")
        s1x_time_path_pdf = os.path.join(plot_dir, f"{s1x_time_base}.pdf")
        plt.savefig(s1x_time_path_png, dpi=600, bbox_inches="tight", pad_inches=0.02)
        plt.savefig(s1x_time_path_pdf, dpi=600, bbox_inches="tight", pad_inches=0.02)
        print(f"Comparison figure (s1x) saved as: {s1x_time_path_png}, {s1x_time_path_pdf}")

        plt.close(fig_s1x_time)

        # 3. spin evolution
        fig_ax1 = plt.figure(figsize=(12, 8))  

        fig_width2, fig_height2 = fig_ax1.get_size_inches()


        ax1_separate = plt.gca()


        line_sx, = ax1_separate.plot(t, sx, color=colors['blue'], linewidth=0.6, alpha=0.8, label='$s_x$')
        line_sy, = ax1_separate.plot(t, sy, color=colors['red'], linewidth=0.6, alpha=0.8, label='$s_y$')
        line_sz, = ax1_separate.plot(t, sz, color=colors['green'], linewidth=0.6, alpha=0.8, label='$s_z$')
        line_s0x_avg, = ax1_separate.plot(t_x_average, sx_average, color='cyan', linestyle='--', linewidth=2.0, label='$s_{0x}$ (avg)')
        line_s0y_avg, = ax1_separate.plot(t_y_average, sy_average, color='magenta', linestyle='--', linewidth=2.0, label='$s_{0y}$ (avg)')
        line_s0z_avg, = ax1_separate.plot(t_z_average, sz_average, color='blue', linestyle='--', linewidth=2.0, label='$s_{0z}$ (avg)')
        line_s0x_th, = ax1_separate.plot(t_trace, sx0_tm, color=colors['orange'], linewidth=2.0, label='$s_{0x}$ (theory)')
        line_s0y_th, = ax1_separate.plot(t_trace, sy0_tm, color=colors['purple'], linewidth=2.0, label='$s_{0y}$ (theory)')
        line_s0z_th, = ax1_separate.plot(t_trace, sz0_tm, color=colors['brown'], linewidth=2.0, label='$s_{0z}$ (theory)')

        ax1_separate.set_xlabel('$t$', fontsize=32, labelpad=8)
        ax1_separate.set_ylabel('Spin Components', fontsize=32, labelpad=8)
        ax1_separate.text(0.02, 0.98, '(a)', transform=ax1_separate.transAxes, fontsize=34, fontweight='bold', 
                        verticalalignment='top', horizontalalignment='left')
        ax1_separate.grid(True, alpha=0.3, linestyle='--')
        ax1_separate.tick_params(axis='both', which='major', width=2.0, length=10, labelsize=28, pad=6)

        ax1_twin = ax1_separate.twinx()
        line_lz, = ax1_twin.plot(t_trace, y_rk4[:,2], color=colors['cyan'], linestyle='--', linewidth=2.0, alpha=0.8, label='$L_z$')
        ax1_twin.set_ylabel('$L_z$', fontsize=32, color=colors['cyan'], labelpad=8)
        ax1_twin.tick_params(axis='y', labelcolor=colors['cyan'], labelsize=28, width=2.0, length=10, pad=6)
        leg = ax1_separate.legend(loc='upper right', fontsize=22, frameon=True, framealpha=0.5,
                                edgecolor='black', fancybox=True, ncol=3,  
                                bbox_to_anchor=(0.82, -0.17),  
                                handlelength=2.0, handletextpad=0.5, columnspacing=0.8)

        plt.subplots_adjust(left=0.05, right=0.75, top=0.95, bottom=0.08)  

        ax1_base = f"{base_filename}_spin_components_evolution"
        ax1_path_png = os.path.join(plot_dir, f"{ax1_base}.png")
        ax1_path_pdf = os.path.join(plot_dir, f"{ax1_base}.pdf")
        plt.savefig(ax1_path_png, dpi=600, bbox_inches="tight", pad_inches=0.02)
        plt.savefig(ax1_path_pdf, dpi=600, bbox_inches="tight", pad_inches=0.02)
        print(f"Separate spin components evolution figure saved as: {ax1_path_png}, {ax1_path_pdf}")

        #plt.show()
