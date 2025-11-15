import glob, os, re, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline
from scipy.integrate import quad
import scienceplots
from classy import Class
import warnings

scales = np.array([0.25, 0.29, 0.33, 0.4, 0.5, 0.67, 1.])
redshifts = 1/scales - 1

class PowerSpectrumAnalyzer:
    def __init__(self, data_dir=None):
        """
        Initialize and automatically read power spectrum data from files.

        Parameters:
        - data_dir: Directory containing the power spectrum files (defaults to current directory).
        """
        self.data_dir = data_dir if data_dir else '/home/vpedre/mochi_CONCEPT/output/leftraru_hi'
        self.k, self.P_sims, self.P_corr, self.P_lins = {}, {}, {}, {}
        self.sigma8_sims, self.sigma8_corr, self.sigma8_lins = {}, {}, {}
        
        self._read_data()

    def a2z(self, a):
        return 1/a-1
    
    def z2a(self, z):
        return 1/(1+z)

    def _read_data(self):
        """Reads power spectrum data from files in the specified directory."""
        file_pattern = f'{self.data_dir}/powerspec_*'

        for filename in sorted(glob.glob(file_pattern), key=os.path.getmtime):
            if filename.endswith('.png'):
                continue

            base = os.path.basename(filename)

            # Match key=value entries (e.g., H0=..., Ωb=...)
            param_matches = re.findall(r'([^_=]+)=([^_]+)', base)
            if not param_matches:
                print(f"Skipping {filename}: no parameters found.")
                continue

            # Parse values and build a key
            parsed_params = {}
            for key, val in param_matches:
                try:
                    parsed_params[key] = float(val)
                except ValueError:
                    parsed_params[key] = val  # fallback to string

            key = tuple(sorted(parsed_params.items()))

            try:
                with open(filename, "r") as f:
                    lines = f.readlines()

                sigma8_values = None
                for line in lines:
                    if re.search(r"σ₈\s*=", line):
                        sigma8_values = list(map(float, re.findall(r"[\d\.e\+\-]+", line)))
                        break

                k_, P_sim, P_c, P_lin = np.loadtxt(filename, usecols=(0, 2, 3, 4), unpack=True)
                mask = ~np.isnan(P_lin)

                self.P_sims[key] = P_sim[mask]
                self.P_corr[key] = P_c[mask]
                self.P_lins[key] = P_lin[mask]
                self.k[key] = k_[mask]

                if sigma8_values:
                    self.sigma8_sims[key] = sigma8_values[0]
                    self.sigma8_corr[key] = sigma8_values[1]
                    self.sigma8_lins[key] = sigma8_values[2]

            except Exception as e:
                print(f"Error processing {filename}:\n{e}")
                continue

        return self.k, self.P_corr, self.P_lins, self.sigma8_corr, self.sigma8_lins


    
    def class_pk(self, zz=redshifts, class_pars=None):
        model = Class()
        try:
            model.set(class_pars)
            model.compute()

            # Define k values to evaluate the power spectrum at
            kk = np.logspace(-4, np.log10(3), 100)  # k values in h/Mpc

            # Initialize a dictionary to store results for each redshift
            Pk_dict = {}

            # Compute power spectrum for all redshifts in z_array
            for z in zz:
                Pk_dict[z] = [model.pk(k, z) for k in kk]  # Store power spectrum for each z

            return kk, Pk_dict
        finally:
            model.empty()

    def plot_Pk(self, models, norm=True, class_curves=None):
        fig, ax = plt.subplots()
        if isinstance(models, tuple):
            models = [models]

        for model in models:
            if norm == 'read':
                normalization = (self.sigma8_lins[model] / self.sigma8_corr[model])**2
            elif norm == 'compute':
                sigma8_lin = self.compute_sigma8(model, use_corrected=False)
                sigma8_corr = self.compute_sigma8(model, use_corrected=True)
                normalization = (sigma8_lin / sigma8_corr)**2
            elif norm == 'bad':
                normalization = self.P_lins[model][0]/self.P_corr[model][0]
            else:
                normalization = 1
                warnings.warn("Normalization set to 1 by default. Make sure this is intended.", UserWarning)

            if class_curves is not None and model in class_curves:
                k_arr, Pk_dict = class_curves[model]
                a = dict(model)['a']
                z = self.a2z(a)
                if z in Pk_dict:
                    ax.loglog(k_arr, Pk_dict[z], label='CLASS')
            w0 = dict(model)['w0']
            wa = dict(model)['wa']
            label_suffix = rf"$w_0 = {w0},\, w_a = {wa}$"
            ax.loglog(self.k[model], self.P_lins[model], '-.', label=label_suffix)
            ax.loglog(self.k[model], normalization * self.P_corr[model], label=label_suffix)

        a_val = dict(models[0])['a']
        ax.set_ylabel(r'$P\,(k,z={:.1f})$ (Mpc$^{{3}}$)'.format(1/a_val - 1))
        ax.set_xlabel(r'$k$ (Mpc$^{-1}$)')
        ax.legend()
        fig.savefig('plot.jpeg')
    
    def interpolate_Pk(self, model, use_corrected=True):
        k_values = self.k[model]
        P_values = self.P_corr[model] if use_corrected else self.P_lins[model]
        return interp1d(k_values, P_values, kind='cubic', bounds_error=False, fill_value="extrapolate")

    
    def window_function(self, k, R=8):
        kR = k * R
        term = np.sin(kR) / kR - np.cos(kR)
        return 3 / kR**2 * term
    
    def compute_sigma8(self, model, use_corrected=True, R=8):
        pk_interp = self.interpolate_Pk(model, use_corrected)

        def integrand(lk):
            k = np.exp(lk)
            return np.exp(3 * lk) * pk_interp(k) * self.window_function(k, R)**2

        integral, _ = quad(integrand, np.log(self.k[model][0]), np.log(self.k[model][-1]), limit=200)

        return np.sqrt(integral / (2 * np.pi**2))
