import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import find_peaks

class PT_plt:
    def __init__(self, h_file, particle_name = '0'):
        self.h_file = h_file # h_file should be a handle to an opened h5 file created by h5py.File()
        self.particle_name = particle_name
        self.axis_dict = {'t':0, 'zeta':1, 'x':2, 'y':3, 'pz':4, 'px':5, 'py':6}
        self.data = self.h_file[self.particle_name]
        
    def get(self, axis):
        return self.data[:, self.axis_dict[axis]]
        
    def get_half_step(self, axis):
        return (self.data[1:, self.axis_dict[axis]]+self.data[:-1, self.axis_dict[axis]])/2
        
    def get_d_half_step(self, axis):
        return self.data[1:, self.axis_dict[axis]]-self.data[:-1, self.axis_dict[axis]]
        
    def save_gamma(self):
        self.gamma = np.sqrt(np.square(self.get('pz')) + np.square(self.get('px')) + np.square(self.get('py')) + 1.)
        
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
    from scipy.integrate import cumtrapz
    r_e = 2.8179403227e-15*1e5
    A = 1.
    gamma0 = 1e5
    x1 = A*gamma0**-0.25
    p1 = x1*(gamma0/2)**.5
    omega = (2*gamma0)**-0.5
    gamma_z0 = gamma0/(1+x1*x1*gamma0/4)**0.5
    gamma_w = gamma_z0
    beta_w = (1-gamma_w**-2)**0.5
    with h5py.File('output.h5','r') as h_file:
        plotter = PT_plt(h_file)
        fig,ax = plt.subplots(figsize=(4, 3))
        plt.plot(plotter.get('t'), plotter.get('zeta')+0.0014851849568716428, 'k-', label='$\\zeta$')
        #plotter.save_gamma()
        #plt.plot(plotter.get_half_step('t'), plotter.get_d_gamma_half_step()/plotter.get_d_half_step('t'), 'r--', label='$\\dot{{\\gamma}}$')
        plt.xlabel('$t$')
        #plt.ylabel('$\\zeta$, $\\dot{{\\gamma}}$')
        plt.ylabel('$\\zeta-\\zeta_{{0l}}$')
        #plt.plot([0.,300.],[8.84e-6,8.84e-6],'k-',[0.,300.],[-8.84e-6,-8.84e-6],'k-',linewidth=0.5)
        #plt.plot([0.,300.],[1.33e-5,1.33e-5],'r-',[0.,300.],[-1.33e-5,-1.33e-5],'r-',linewidth=0.5)
        plt.plot([0.,1000.],[8.84e-7,8.84e-7],'k-',[0.,1000.],[-8.84e-7,-8.84e-7],'k-',linewidth=0.5)
        '''plt.plot(plotter.get('t'), plotter.get('x'), 'k-', label='$\\zeta$')
        plt.plot(plotter.get('t'), (plotter.save_get_gamma())**-.25, 'r--', label='$A\\gamma^{{-1/4}}$')'''
        '''zeta=plotter.get('zeta')
        x=plotter.get('x')
        fig,ax = plt.subplots()
        plt.plot(zeta, x, 'k-')
        plt.arrow(zeta[31],x[31],zeta[34]-zeta[28],x[34]-x[28],width=2e-8,head_length=3e-3,head_width=1.5e-7,length_includes_head=True,overhang=0.5,color='r')
        plt.arrow(zeta[164],x[164],zeta[167]-zeta[161],x[167]-x[161],width=2e-8,head_length=3e-3,head_width=1.5e-7,length_includes_head=True,overhang=0.5,color='r')
        plt.arrow(zeta[230],x[230],zeta[233]-zeta[227],x[233]-x[227],width=2e-8,head_length=3e-3,head_width=1.5e-7,length_includes_head=True,overhang=0.5,color='r')
        plt.arrow(zeta[370],x[370],zeta[373]-zeta[367],x[373]-x[367],width=2e-8,head_length=3e-3,head_width=1.5e-7,length_includes_head=True,overhang=0.5,color='r')
        #plt.arrow(zeta[0],x[0],zeta[4042]-zeta[0],x[4042]-x[0],width=2e-7,head_length=3e-7,head_width=1.5e-8,length_includes_head=True,overhang=0.5,color='b')
        ax.annotate("", xy=(zeta[0],x[0]), xytext=(zeta[4042],x[4042]), arrowprops=dict(arrowstyle="->", head_width=0.2))
        #plt.scatter(zeta[4045:],x[4045:])
        plt.xlabel('$k_p\\zeta$')
        plt.ylabel('$k_p x$')'''
        '''gamma = plotter.save_get_gamma()
        peak_ind = plotter.find_peak_inds('x')
        t = plotter.get('t')
        #plt.plot(t[peak_ind],np.abs(plotter.get('x')[peak_ind]),'k-')
        #plt.plot(t,A*gamma**-.25,'r--')
        x = plotter.get('x')
        px = plotter.get('px')
        l1=403
        n2=int(len(x)/3.44)
        l2=556
        l3=801
        fig,ax = plt.subplots(figsize=(4, 3))
        #plt.figure()
        plt.plot(x[:l1],px[:l1],'k',label='$t_0$')#'$\\omega_p t=0$')
        plt.plot(x[n2:n2+l2],px[n2:n2+l2],'r',label='$t_1$')#'$\\omega_p t= 10^5$')
        plt.plot(x[-l3:],px[-l3:],'b',label='$t_2$')#'$\\omega_p t=3.5\\times 10^5$')
        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xlabel('$k_p x$', loc="right")
        plt.ylabel('$p_x/m_ec$', loc="top", rotation=0)
        plt.xlim([None,0.17])
        plt.ylim([None,8.])
        print(t[l1//2])
        print(t[n2+l2//2])
        print(t[-l3//2])'''
        #plotter.save_gamma()
        '''plt.plot(plotter.get('t'), plotter.get('zeta'), 'k-', label='$\\zeta$')
        plt.plot(plotter.get('t'), plotter.gamma, 'r-', label='$\\gamma$')
        plt.xlabel('$t$')
        plt.ylabel('$\\gamma$')
        plt.legend()
        plt.figure(figsize=(4, 3))'''

        '''n_dot_plot = 128
        S, t = plotter.get_S()
        inds = np.linspace(0, len(t)-1, n_dot_plot ,dtype=int)
        t = t[inds]
        S = S[inds]
        plt.plot(t, S, 'k-', linewidth=.5, label='Numerical')
        Dt = t - t[0]
        plt.plot(t, 1/(r_e/8/2**0.5/np.pi*gamma0**0.5*Dt + 1/S[0]), 'r--', label=' ')
        cumint_sqrt_gamma = cumtrapz(np.sqrt(plotter.gamma), plotter.get('t'))
        t = plotter.get_half_step('t')
        S = 1./(r_e/8/2**0.5/np.pi*cumint_sqrt_gamma + 1/S[0])
        inds = np.linspace(0, len(t)-1, n_dot_plot ,dtype=int)
        plt.plot(t[inds], S[inds], 'b--', label=' ')'''
        '''#plt.plot(plotter.get('t'), plotter.get('pz'), 'k-', label='$p_z$')
        #plt.plot(plotter.get('t'), plotter.gamma, 'k-', label='$\\gamma$')
        #plt.plot(plotter.get('t'), plotter.get_gamma()-plotter.get('pz'), 'k-', label='$\\gamma-p_z$')
        #plt.plot(plotter.get_half_step('t'), (plotter.get_d_half_step('zeta')+beta_w)*plotter.get_gamma_half_step()-plotter.get_half_step('pz'), 'r-', label='$p_z$')
        plt.plot(plotter.get_half_step('t'), plotter.get_d_half_step('pz')+.5*plotter.get_half_step('zeta'), 'r-', label='$\\dot{p}_z$')'''
        #plt.legend(loc=4)
        plt.minorticks_on()
        plt.tight_layout()
        plt.show()
