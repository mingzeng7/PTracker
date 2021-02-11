import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import find_peaks

class PT_plt:
    def __init__(self, h_file, particle_name = '0'):
        self.h_file = h_file # h_file should be a handle to an opened h5 file created by h5py.File()
        self.particle_name = particle_name
        self.axis_dict = {'t':0, 'zeta':1, 'x':2, 'pz':3, 'px':4}
        self.data = self.h_file[self.particle_name]
        
    def get(self, axis):
        return self.data[:, self.axis_dict[axis]]
        
    def get_half_step(self, axis):
        return (self.data[1:, self.axis_dict[axis]]+self.data[:-1, self.axis_dict[axis]])/2
        
    def get_d_half_step(self, axis):
        return self.data[1:, self.axis_dict[axis]]-self.data[:-1, self.axis_dict[axis]]
        
    def save_gamma(self):
        self.gamma = np.sqrt(np.square(self.get('pz')) + np.square(self.get('px')) + 1.)
        
    def save_get_gamma(self):
        self.save_gamma()
        return self.gamma
        
    def get_gamma_half_step(self):
        return (self.gamma[1:]+self.gamma[:-1])/2
        
    def get_d_gamma_half_step(self):
        return self.gamma[1:]-self.gamma[:-1]
        
    def find_peak_inds(self, axis):
        return find_peaks(np.abs(self.get(axis)))[0]
        
    def get_x1_times_px1(self):
        t = self.get('t')
        x_peak_inds = self.find_peak_inds('x')
        x_peaks = np.abs(self.get('x')[x_peak_inds])
        px_peak_inds = self.find_peak_inds('px')
        px_peaks = np.abs(self.get('px')[px_peak_inds])
        t_x_peaks = t[x_peak_inds]
        px_peaks_interp = np.interp(t_x_peaks, t[px_peak_inds], px_peaks)
        return x_peaks*px_peaks_interp, t_x_peaks
        
    def get_S(self):
        x1_times_px1, t = self.get_x1_times_px1()
        return x1_times_px1*np.pi, t
        
if __name__ == "__main__":
    r_e = 2.8179403227e-10
    A = 1.
    gamma0 = 1e5
    x1 = A*gamma0**-0.25
    p1 = x1*(gamma0/2)**.5
    omega = (2*gamma0)**-0.5
    with h5py.File('output.h5','r') as h_file:
        plotter = PT_plt(h_file)
        plt.plot(plotter.get('t'), plotter.get('zeta'), 'k-', label='$\\zeta$')
        plotter.save_gamma()
        plt.plot(plotter.get('t'), plotter.gamma, 'r-', label='$\\gamma$')
        plt.figure(figsize=(4, 3))

        S, t = plotter.get_S()
        plt.plot(t, S, 'k-', linewidth=.5, label='Numerical')
        Dt = t - t[0]
        plt.plot(t, 1/(r_e/8/2**0.5/np.pi*gamma0**0.5*Dt + 1/S[0]), 'r--', label=' ')
        plt.ylabel('$S$')
        plt.xlabel('$t$')
        plt.minorticks_on()
        plt.legend(loc=1)
        plt.tight_layout()
        plt.show()
