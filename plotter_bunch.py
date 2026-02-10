import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import traceback
import pylibconfig2 as libconfig
from tqdm import tqdm
from matplotlib.animation import FuncAnimation, FFMpegWriter
import matplotlib.animation as animation
from datetime import datetime

# Create Config object
config = libconfig.Config()

# Read configuration file
config.read_file('input.cfg')

class BunchAnalyzer:
    def __init__(self, filename, chunk_size=1000, max_time=None):
        self.filename = filename
        self.time = None
        self.attrs = {}
        self.var_names = []
        self.var_indices = {}
        self.chunk_size = chunk_size  
        self.max_time = max_time  
        self.load_config_parameters()
        self.load_metadata()         
        self._cache = {}
    
    def load_config_parameters(self):
        try:
            wake_config = config.lookup('wake')
            self.kappa_square = wake_config.get('kappa_square', 0.5)
            bunch_config = config.lookup('particle.bunch')
            self.gamma0_mean = bunch_config.get('gamma0_mean', 10.0)
            self.if_RR1 = config.lookup('if_RR1')
            self.if_SG = config.lookup('if_SG')
            self.if_SG_grad = config.lookup('if_SG_grad')

            self.elec_anomaly = 1.1596521884e-3

            self.sampling_factor = config.lookup('sampling_factor')
            self.t1 = config.lookup('t1')
            self.particle = config.lookup('particle')
            self.kappa_square = 0.5
            wake = config.lookup('wake')
            self.kp_SI = wake._lookup(['kp_SI'])
            self.gama_w = wake._lookup(['gamma_w'])
            self.lambd = wake._lookup(['lambda'])
            self.f_z0 = wake._lookup(['f_z0'])
            print(f"Parameters from .cfg file:")
            print(f"  kappa_square = {self.kappa_square}")
            print(f"  gamma0_mean = {self.gamma0_mean}")
            
        except Exception as e:
            print(f"Failed to load .cfg file: {e}")
            # Set default values
            self.kappa_square = 0.5
            self.gamma0_mean = 10.0
            print("Using default parameters")
    
    def load_metadata(self):
        print(f"Loading file metadata: {self.filename}")
        try:
            with h5py.File(self.filename, 'r') as f:
                # Read time data - complete time series
                self.time = f['time'][:]
                
                # If max analysis time is set, find corresponding timestep indices
                if self.max_time is not None:
                    self.analysis_timesteps = np.where(self.time <= self.max_time)[0]
                    if len(self.analysis_timesteps) == 0:
                        print(f"Warning: Max time {self.max_time} is smaller than first timestep, using all data")
                        self.analysis_timesteps = np.arange(len(self.time))
                    else:
                        self.analysis_timesteps = np.arange(self.analysis_timesteps[-1] + 1)
                    print(f"Analysis time range: 0 to {self.max_time}, timesteps: {len(self.analysis_timesteps)}")
                else:
                    self.analysis_timesteps = np.arange(len(self.time))
                    print(f"Analyzing full time range, timesteps: {len(self.analysis_timesteps)}")
                
                # Read dataset shape information
                dataset = f['bunch_data']
                self.total_particles = dataset.shape[0]
                self.total_timesteps = dataset.shape[1]
                self.total_variables = dataset.shape[2]
                
                print(f"Data shape: particles={self.total_particles}, timesteps={self.total_timesteps}, variables={self.total_variables}")
                print(f"Simulation duration: {self.time[-1]:.2f}")
                
                # Read attributes
                print("Dataset attributes:")
                for attr_name in dataset.attrs:
                    attr_value = dataset.attrs[attr_name]
                    if isinstance(attr_value, bytes):
                        self.attrs[attr_name] = attr_value.decode('utf-8')
                    else:
                        self.attrs[attr_name] = attr_value
                    print(f"  {attr_name}: {self.attrs[attr_name]}")
                
                # Build variable name mapping
                self.var_names = ['t', 'zeta', 'x', 'y', 'pz', 'px', 'py', 'sz', 'sx', 'sy']
                self.var_indices = {name: idx for idx, name in enumerate(self.var_names)}
                print(f"Variable mapping: {self.var_indices}")
                
        except Exception as e:
            print(f"Failed to load metadata: {e}")
            traceback.print_exc()
            raise
    
    def get_timestep_data(self, timestep_idx):
        """Get data for a single timestep"""
        with h5py.File(self.filename, 'r') as f:
            dataset = f['bunch_data']
            return dataset[:, timestep_idx, :]
    
    def get_timestep_range_data(self, start_idx, end_idx):
        """Get data for a range of timesteps"""
        with h5py.File(self.filename, 'r') as f:
            dataset = f['bunch_data']
            return dataset[:, start_idx:end_idx, :]
    
    def calculate_theoretical_polarization_direct(self):
        """
        Calculate theoretical polarization change
        """
        # Electron anomalous magnetic moment
        a_e = 0.00115965218128
        
        print("Calculating theoretical polarization...")
        
        # Get initial data
        initial_data = self.get_timestep_data(0)
        idx = self.var_indices
        
        # Get initial polarization components
        Px_0 = np.mean(initial_data[:, idx['sx']])
        Py_0 = np.mean(initial_data[:, idx['sy']])
        Pz_0 = np.mean(initial_data[:, idx['sz']])
        
        print(f"Initial polarization components: Px_0={Px_0:.6f}, Py_0={Py_0:.6f}, Pz_0={Pz_0:.6f}")
        
        # Calculate initial physical quantities
        pz_initial = initial_data[:, idx['pz']]
        px_initial = initial_data[:, idx['px']]
        py_initial = initial_data[:, idx['py']]
        gamma0_i = np.sqrt(1 + pz_initial**2 + px_initial**2 + py_initial**2)
        
        # Calculate initial angular momentum
        x_initial = initial_data[:, idx['x']]
        y_initial = initial_data[:, idx['y']]
        l_z_initial = x_initial * py_initial - y_initial * px_initial
        
        # Initialize accumulated phase for each particle
        num_particles = len(gamma0_i)
        accumulated_phase = np.zeros(num_particles)
        
        # Calculate ω_T for initial timestep
        valid_indices = gamma0_i > 0
        prev_omega_T = np.zeros(num_particles)
        prev_omega_T[valid_indices] = (
            0.5 * (gamma0_i[valid_indices] * self.kappa_square + self.f_z0**2) * 
            (a_e + 1.0/gamma0_i[valid_indices])**2 - 
            self.kappa_square * (a_e + 1.0/gamma0_i[valid_indices])
        ) / gamma0_i[valid_indices] * l_z_initial[valid_indices]
        
        # Initialize result array
        num_timesteps = len(self.analysis_timesteps)
        theoretical_polarization = np.zeros(num_timesteps)
        time_array = self.time[self.analysis_timesteps]
        
        # Process first timestep
        theoretical_polarization[0] = np.sqrt(Px_0**2 + Py_0**2 + Pz_0**2)
        
        # Process subsequent timesteps
        for i in tqdm(range(1, num_timesteps), desc="Calculating theoretical polarization"):
            # Get current timestep data
            timestep_idx = self.analysis_timesteps[i]
            current_data = self.get_timestep_data(timestep_idx)
            
            # Calculate particle energy for current timestep
            pz_current = current_data[:, idx['pz']]
            px_current = current_data[:, idx['px']]
            py_current = current_data[:, idx['py']]
            gamma_i = np.sqrt(1 + pz_current**2 + px_current**2 + py_current**2)
            
            # Calculate angular momentum for current timestep
            x_current = current_data[:, idx['x']]
            y_current = current_data[:, idx['y']]
            l_z_current = x_current * py_current - y_current * px_current

            # Calculate ω_T for current timestep (vectorized)
            valid_indices = gamma_i > 0
            current_omega_T = np.zeros(num_particles)
            current_omega_T[valid_indices] = (
                0.5 * (gamma_i[valid_indices] * self.kappa_square + self.f_z0**2) * 
                (a_e + 1.0/gamma_i[valid_indices])**2 - 
                self.kappa_square * (a_e + 1.0/gamma_i[valid_indices])
            ) / gamma_i[valid_indices] * l_z_current[valid_indices]
            
            # Calculate timestep
            dt = time_array[i] - time_array[i-1]
            
            # Update accumulated phase using trapezoidal rule
            phase_increment = 0.5 * (prev_omega_T + current_omega_T) * dt
            accumulated_phase += phase_increment
            
            # Update ω_T from previous timestep
            prev_omega_T = current_omega_T
            
            # Calculate averages of cos(phase) and sin(phase)
            cos_terms = np.cos(accumulated_phase)
            sin_terms = np.sin(accumulated_phase)
            
            # Calculate polarization
            Px_current = np.mean(Px_0*cos_terms - Py_0*sin_terms)
            Py_current = np.mean(Px_0*sin_terms + Py_0*cos_terms)
            
            theoretical_polarization[i] = np.sqrt(Px_current**2 + Py_current**2 + Pz_0**2)
        
        return theoretical_polarization
    
    def calculate_beam_parameters_batch(self):
        """Calculate beam parameters in batches - full time resolution"""
        # Check cache
        if 'beam_params' in self._cache:
            return self._cache['beam_params']
            
        print("Calculating beam parameters in batches...")
        
        num_timesteps = len(self.analysis_timesteps)
        idx = self.var_indices
        
        # Initialize result arrays - full time resolution
        mean_energy = np.zeros(num_timesteps)
        mean_sx = np.zeros(num_timesteps)
        mean_sy = np.zeros(num_timesteps)
        mean_sz = np.zeros(num_timesteps)
        polarization = np.zeros(num_timesteps)
        
        # Process in batches - full time series
        for chunk_start in tqdm(range(0, num_timesteps, self.chunk_size), desc="Processing timesteps"):
            chunk_end = min(chunk_start + self.chunk_size, num_timesteps)
            chunk_length = chunk_end - chunk_start
            
            # Calculate actual timestep indices
            actual_start = self.analysis_timesteps[chunk_start]
            actual_end = self.analysis_timesteps[chunk_end-1] + 1 if chunk_end < num_timesteps else self.analysis_timesteps[-1] + 1
            
            # Read current batch data
            chunk_data = self.get_timestep_range_data(actual_start, actual_end)
            
            # Calculate gamma
            gamma_chunk = np.sqrt(1 + chunk_data[:, :, idx['pz']]**2 + 
                                chunk_data[:, :, idx['px']]**2 + 
                                chunk_data[:, :, idx['py']]**2)
            
            # Calculate energy parameters
            mean_energy_chunk = np.mean(gamma_chunk, axis=0)
            
            # Calculate average spin
            mean_sx_chunk = np.mean(chunk_data[:, :, idx['sx']], axis=0)
            mean_sy_chunk = np.mean(chunk_data[:, :, idx['sy']], axis=0)
            mean_sz_chunk = np.mean(chunk_data[:, :, idx['sz']], axis=0)
            
            # Calculate spin polarization
            polarization_chunk = np.sqrt(
                mean_sx_chunk**2 + 
                mean_sy_chunk**2 + 
                mean_sz_chunk**2
            )
            
            # Ensure shape matches
            if len(mean_energy_chunk) == chunk_length:
                mean_energy[chunk_start:chunk_end] = mean_energy_chunk
                mean_sx[chunk_start:chunk_end] = mean_sx_chunk
                mean_sy[chunk_start:chunk_end] = mean_sy_chunk
                mean_sz[chunk_start:chunk_end] = mean_sz_chunk
                polarization[chunk_start:chunk_end] = polarization_chunk
            else:
                print(f"Warning: Data shape mismatch - mean_energy_chunk: {mean_energy_chunk.shape}, expected length: {chunk_length}")
                # If shape doesn't match, process each timestep with loop
                for i in range(min(chunk_length, len(mean_energy_chunk))):
                    actual_idx = chunk_start + i
                    if actual_idx < num_timesteps:
                        mean_energy[actual_idx] = mean_energy_chunk[i]
                        mean_sx[actual_idx] = mean_sx_chunk[i]
                        mean_sy[actual_idx] = mean_sy_chunk[i]
                        mean_sz[actual_idx] = mean_sz_chunk[i]
                        polarization[actual_idx] = polarization_chunk[i]
        
        beam_params = {
            'time': self.time[self.analysis_timesteps],
            'mean_energy': mean_energy,
            'mean_sx': mean_sx,
            'mean_sy': mean_sy,
            'mean_sz': mean_sz,
            'polarization': polarization
        }
        
        # Cache results
        self._cache['beam_params'] = beam_params
        return beam_params
    
    def get_all_analysis_data(self):
        """Get all analysis data to avoid repeated calculations"""
        if 'all_data' in self._cache:
            return self._cache['all_data']
            
        print("Starting calculation of all analysis data...")
        
        # Calculate all required parameters
        beam_params = self.calculate_beam_parameters_batch()
        
        all_data = {
            'beam_params': beam_params,
        }
        
        self._cache['all_data'] = all_data
        return all_data
    
    def calculate_theoretical_polarization_from_data(self):
        """Calculate theoretical polarization from existing data"""
        return self.calculate_theoretical_polarization_direct()
    
    def plot_polarization_comparison(self, save_path=None, precomputed_theoretical_polarization=None):
        """Plot polarization evolution comparison"""
        # Get data
        data = self.get_all_analysis_data()
        beam_params = data['beam_params']
        
        # Use precomputed theoretical polarization or recalculate
        if precomputed_theoretical_polarization is not None:
            theoretical_polarization = precomputed_theoretical_polarization
        else:
            theoretical_polarization = self.calculate_theoretical_polarization_from_data()
        
        print("Plotting polarization comparison...")
        
        # Set large font style
        plt.rcParams.update({
            'font.size': 40,
            'font.family': 'serif',
            'mathtext.fontset': 'stix',
            'axes.labelsize': 40,
            'axes.titlesize': 40,
            'legend.fontsize': 40,
            'xtick.labelsize': 40,
            'ytick.labelsize': 40
        })
        
        # Create figure - using square size
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot polarization evolution comparison
        line1 = ax.plot(beam_params['time'], beam_params['polarization'], 
                color='#1f77b4', linewidth=3, label='P (numerical)')[0]
        line2 = ax.plot(beam_params['time'], theoretical_polarization, 
                'red', linewidth=3, linestyle='--', 
                label='P (theory)')[0]
        
        # Create twin axis and plot gamma
        ax_twin = ax.twinx()
        line3 = ax_twin.plot(beam_params['time'], beam_params['mean_energy'], 
                color='green', linewidth=3, label=r'$\gamma$')[0]
        
        # Set axis labels
        ax.set_xlabel('t', fontsize=40)
        ax.set_ylabel('P', fontsize=40)
        ax_twin.set_ylabel(r'$\gamma$', fontsize=40)
        
        # Set x-axis to scientific notation
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
        
        # Set gamma axis to scientific notation
        ax_twin.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        # Set scientific notation font size
        ax.xaxis.get_offset_text().set_fontsize(40)
        if ax_twin.yaxis.get_offset_text():
            ax_twin.yaxis.get_offset_text().set_fontsize(40)
                
        # Merge legends - put all three lines in same legend
        lines = [line1, line2, line3]
        labels = [line.get_label() for line in lines]
        
        # Create merged legend - place inside plot at upper right
        ax.legend(lines, labels, loc='lower right', fontsize=40)
        
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Polarization comparison plot saved to: {save_path}")
        
        plt.show()
    
    def generate_report(self, precomputed_theoretical_polarization=None):
        """Generate analysis report"""
        # Get data
        data = self.get_all_analysis_data()
        beam_params = data['beam_params']
        
        # Use precomputed theoretical polarization or recalculate
        if precomputed_theoretical_polarization is not None:
            theoretical_polarization = precomputed_theoretical_polarization
        else:
            theoretical_polarization = self.calculate_theoretical_polarization_from_data()
        
        print("="*70)
        print(f"Bunch Dynamics Analysis Report (Direct Method for Theoretical Polarization)")
        print("="*70)
        print(f"Data file: {self.filename}")
        print(f"Particle count: {self.total_particles}")
        print(f"Timestep count: {self.total_timesteps}")
        print(f"Simulation duration: {self.time[-1]:.2f}")
        if self.max_time is not None:
            print(f"Analysis time range: 0 to {self.max_time}")
            print(f"Analysis timesteps: {len(self.analysis_timesteps)}")
        
        print("\nTheoretical parameters:")
        print(f"  kappa_square: {self.kappa_square}")
        print(f"  gamma0_mean: {self.gamma0_mean}")
        print(f"  a_e (electron anomalous magnetic moment): 0.00115965218128")
        print(f"  f_z0: {self.f_z0}")
        print(f"  Initial polarization: {beam_params['polarization'][0]:.4f}")
        
        print("\nEnergy evolution:")
        print(f"  Initial average energy: {beam_params['mean_energy'][0]:.2f}")
        print(f"  Final average energy: {beam_params['mean_energy'][-1]:.2f}")
        
        print("\nSpin polarization:")
        print(f"  Initial polarization (numerical): {beam_params['polarization'][0]:.6f}")
        print(f"  Final polarization (numerical): {beam_params['polarization'][-1]:.6f}")
        print(f"  Final polarization (theoretical): {theoretical_polarization[-1]:.6f}")
        print(f"  Polarization difference (absolute): {abs(beam_params['polarization'][-1] - theoretical_polarization[-1]):.6f}")
        
        print("\nAnalysis completed!")


def IF_RR(if_RR):
    if if_RR:
        return 1.
    else:
        return 0 

def create_plot_output_path(output_file, plot_type, timestamp, file_format="pdf"):
    """
    Create plot output path
    
    Parameters:
    output_file: Original data file path
    plot_type: Plot type, e.g., "polarization_comparison"
    timestamp: Timestamp
    file_format: Output file format, default "pdf"
    
    Returns:
    New output path
    """
    import os
    
    # Get original file directory and base filename
    original_dir = os.path.dirname(output_file)
    original_basename = os.path.basename(output_file)
    original_name_no_ext = os.path.splitext(original_basename)[0]
    
    # Create new output directory
    if original_dir:
        # If original directory contains path separator, may need to extract directory name
        # e.g., "data_20251228" -> "data_20251228_plot"
        new_dir = f"{original_dir}_plot"
    else:
        new_dir = "plots"
    
    # Ensure new directory exists
    os.makedirs(new_dir, exist_ok=True)
    
    # Build new output filename
    new_filename = f"{original_name_no_ext}_{plot_type}.{file_format}"
    
    # Return full path
    return os.path.join(new_dir, new_filename)

def main():
    output_file = config.lookup('output_file')
    filename = output_file
    if_RR = config.lookup('if_RR')
    r_e = 2.8179403262e-15 * 1e5 * IF_RR(if_RR)  # 1e5 is kp_si

    print("Starting analysis...")
    print(f"Current working directory: {os.getcwd()}")
    
    # Check if file exists
    if not os.path.exists(filename):
        print(f"\nError: File '{filename}' does not exist!")
        return
    
    # Backup input.cfg file content as txt
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    try:
        with open('input.cfg', 'r') as cfg_file:
            cfg_content = cfg_file.read()
        
        # Create backup filename using same format as images
        cfg_backup_filename = f"{output_file.replace('.h5', '')}_config_backup_{timestamp}.txt"
        
        # Write backup file
        with open(cfg_backup_filename, 'w') as backup_file:
            backup_file.write(f"Configuration Backup - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            backup_file.write("=" * 50 + "\n\n")
            backup_file.write(cfg_content)
        
        print(f"Configuration backup saved as: {cfg_backup_filename}")
    except Exception as e:
        print(f"Warning: Could not backup configuration file. Error: {e}")
    
    try:
        # Use full time resolution version
        print(f"\nCreating full time resolution analyzer...")
        
        # Option to limit analysis time range
        max_time = None  # e.g., only analyze first 5e4 time units
        
        if max_time is not None:
            print(f"Limiting analysis time range to 0 to {max_time}")
        
        analyzer = BunchAnalyzer(
            filename, 
            chunk_size=1000,      # Batch processing size
            max_time=max_time     # Maximum analysis time
        )
        
        print("Calculating theoretical polarization using direct method...")
        
        # Precalculate all analysis data
        print("Precalculating all analysis data...")
        all_data = analyzer.get_all_analysis_data()
        
        # Calculate theoretical polarization
        print("Calculating theoretical polarization...")
        theoretical_polarization = analyzer.calculate_theoretical_polarization_from_data()
        
        # Plot polarization comparison
        print(f"\nPlotting polarization comparison (large font)...")
        plot_save_path = create_plot_output_path(
            output_file=output_file,
            plot_type="polarization_comparison_direct_method",
            timestamp=timestamp,
            file_format="pdf"
        )
        analyzer.plot_polarization_comparison(save_path=plot_save_path, precomputed_theoretical_polarization=theoretical_polarization)
        
        # Generate analysis report
        print(f"\nGenerating analysis report...")
        analyzer.generate_report(precomputed_theoretical_polarization=theoretical_polarization)
        
        print("Analysis completed!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()