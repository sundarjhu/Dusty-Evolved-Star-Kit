import os,math,subprocess,pdb
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from fnmatch import fnmatch
from multiprocessing import cpu_count, Pool
from functools import partial
from astropy.table import Table, Column
from matplotlib import rc
from scipy import interpolate

rc('text', usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

# options
directory = 'year4'
distance_in_kpc = 8
assumed_rgd = 200.0000

# set variables
targets = []
target_names = []

# solar constant = 1379 W
# distance to sun in kpd 4.8483E-9
distance_norm = math.log10(((int(distance_in_kpc)/4.8482E-9)**2)/1379)
min_norm = 1e-14
max_norm = 1e-10
ntrials = 1000
trials = np.linspace(min_norm, max_norm, ntrials)
latex_array = []

# for multiple sources
for item in os.listdir('./visir_spectra/'):
    if fnmatch(item, "*flux_calibrated.txt"):
        targets.append('visir_spectra/'+item)

# for select sources
targets = ['visir_spectra/IRAS-17030-3053_flux_calibrated.txt']

grid_dusty = Table.read('dusty_models.fits')
grid_outputs = Table.read('dusty_outputs.fits')


def get_data(filename):
    x, y = np.loadtxt(filename, delimiter=',', unpack=1, skiprows=110)
    y = y * u.Jy
    y = y.to(u.W/(u.m * u.m), equivalencies=u.spectral_density(x * u.um))
    return x, np.array(y)


def least2(data, model):
    return np.square(model - data).sum()


def trim(data, model):
    # gets dusty model range that matches data
    indexes = np.where(np.logical_and(model[0] >= np.min(data[0]), model[0] <= np.max(data[0])))
    return np.vstack([model[0][indexes], model[1][indexes]])


def fit_norm(data, model):
    stat_array = []
    # list(trimmed_model[1]*t for t in trials)
    for t in trials:
        stat_values = least2(data[1], trimmed_model[1]*t)
        stat_array.append(stat_values)
    return stat_array


def interpolate_data(model_x,spectra):
    data_wavelength = spectra[0]
    data_flux = spectra[1]
    f = interpolate.interp1d(data_wavelength, data_flux)
    indexes = np.where(np.logical_and(model_x >= np.min(data_wavelength), model_x <= np.max(data_wavelength)))
    resampled_data_flux = f(model_x[indexes])
    return model_x[indexes], resampled_data_flux


# plotting stuff
if len(targets) < 4:
    fig, ax1 = plt.subplots(len(targets), 1, sharex=True, sharey=True, figsize=(16, 20))
else:
    fig, axs = plt.subplots(math.ceil(len(targets)/4), 4, sharex=True, sharey=True, figsize=(16, 20))
    axs = axs.ravel()

for counter, target in enumerate(targets):
    stat_values = []
    raw_data = get_data(target)
    model_x = grid_dusty[0][0][:]
    resampled_data = interpolate_data(model_x, raw_data)
    for model in np.array(grid_dusty):
        trimmed_model = trim(resampled_data, model)
        stat_values.append(fit_norm(resampled_data, model))
    stat_array = np.vstack(stat_values)
    argmin = np.argmin(stat_array)
    model_index = argmin // stat_array.shape[1]
    trial_index = argmin % stat_array.shape[1]

    # calculated luminosity and scales outputs
    luminosity = np.power(10.0, distance_norm - math.log10(trials[trial_index]) * -1)
    scaled_vexp = float(grid_outputs[model_index]['vexp']) * (luminosity / 10000) ** (0.25)
    scaled_mdot = grid_outputs[model_index]['mdot'] * ((luminosity / 10000) ** 0.75) * (assumed_rgd / 200) ** (0.5)

    target_name=(target.split('/')[-1][:15]).replace('IRAS-', 'IRAS ')
    print(target_name)
    latex_array.append((target_name, str(int(round(luminosity/1000))), str(np.round(scaled_vexp, 1)), str(
        int(grid_outputs[model_index]['teff'])), str(int(grid_outputs[model_index]['tinner'])), str(
        grid_outputs[model_index]['odep']), "%.3E" % float(scaled_mdot)))

    # printed output
    print()
    print()
    print(('    Target: '+target_name))
    print()
    print()
    print(("Luminosity\t\t"+str(round(luminosity))))
    print(("Optical depth\t\t"+str(grid_outputs[model_index]['odep'])))
    print(("Expansion velocity (scaled)\t"+str(round(scaled_vexp, 2))))
    print(("Mass loss (scaled)\t\t"+str("%.2E" % float(scaled_mdot))))
    print()
    print()

    # gets data for plotting
    scale_distance = lambda model_val: np.vstack([model_val[0], model_val[1]*trials[trial_index]])
    x_model, y_model = scale_distance(grid_dusty[model_index])
    x_data, y_data = get_data(target)

    # logscale
    x_model = np.log10(x_model)
    y_model = np.log10(y_model)
    x_data = np.log10(x_data)
    y_data = np.log10(y_data)

    # plotting
    if len(targets) == 1:
        ax1.set_xlim(0.21, 1.79)
        ax1.set_ylim(-14.2, -10.51)
        ax1.plot(x_model, y_model, c='k', linewidth=0.7, linestyle='--', zorder=2)
        ax1.plot(x_data, y_data, c='blue')
        ax1.annotate(target_name.replace('-', r'\textendash'), (1.15, -13.75), xycoords='data', fontsize=14)
        ax1.get_xaxis().set_tick_params(which='both', direction='in', labelsize=15)
        ax1.get_yaxis().set_tick_params(which='both', direction='in', labelsize=15)
    else:
        axs[counter].set_xlim(0.21, 1.79)
        axs[counter].set_ylim(-14.2, -10.51)
        axs[counter].plot(x_model, y_model, c='k', linewidth=0.7, linestyle='--', zorder=2)
        axs[counter].plot(x_data, y_data, c='blue')
        axs[counter].annotate(target_name.replace('-', r'\textendash'), (1.15, -13.75), xycoords='data', fontsize=14)
        axs[counter].get_xaxis().set_tick_params(which='both', direction='in', labelsize=15)
        axs[counter].get_yaxis().set_tick_params(which='both', direction='in', labelsize=15)

plt.subplots_adjust(wspace=0, hspace=0)
fig.text(0.5, 0.085, 'log $\lambda$ ($\mu m$)', ha='center', fontsize=16)
fig.text(0.08, 0.5, "log $\lambda$ F$_{\lambda}$ "+"(W m$^{-2}$)", va='center', rotation='vertical', fontsize=16)
fig.savefig('output_seds.png', dpi=500, bbox_inches='tight')
    

a = Table(np.array(latex_array), names=(
    'source', 'L', 'vexp_predicted', 'teff', 'tinner', 'odep', 'mdot'), dtype=(
    'S16', 'int32', 'f8', 'int32', 'int32', 'f8', 'f8'))
a['L'] = a['L']*1000
a.write('output_latex_table.csv', format='csv', overwrite=True)


