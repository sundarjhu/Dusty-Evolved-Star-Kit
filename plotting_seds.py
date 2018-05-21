import os
import math
import subprocess
import pdb
import shutil
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from glob import glob
from fnmatch import fnmatch
from multiprocessing import cpu_count, Pool
from functools import partial
from astropy.table import Table, Column
from matplotlib import rc
from scipy import interpolate
'''
Steve Goldman
Space Telescope Science Institute
May 17, 2018
sgoldman@stsci.edu

This script is for plotting the outputs of the sed_fitting script.
'''

rc('text', usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

input_file = Table.read('fitting_plotting_outputs.csv')
grid_dusty = Table.read('models/'+str(input_file['grid_name'][0])+'_models.fits')
grid_outputs = Table.read('models/'+str(input_file['grid_name'][0])+'_outputs.fits')


def get_data(filename):
    table = Table.read(filename)
    x = np.array(table.columns[0])
    y = np.array(table.columns[1])
    y = y * u.Jy
    y = y.to(u.W/(u.m * u.m), equivalencies=u.spectral_density(x * u.um))
    return x, np.array(y)


def get_supp_data(file):
    supp_array = Table(np.genfromtxt(file, delimiter=',', names=True))
    xs = np.array(supp_array.columns[0])
    ys = np.array(supp_array.columns[1])
    ys = ys * u.Jy
    ys = ys.to(u.W/(u.m * u.m), equivalencies=u.spectral_density(xs * u.um))
    out_supp = np.array(list(zip(xs, ys.value)))
    return out_supp


# plotting stuff
if len(input_file) == 1:
    fig, ax1 = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(16, 20))
elif len(input_file) == 2:
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(16, 20))
elif len(input_file) == 3:
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(16, 20))
else:
    fig, axs = plt.subplots(math.ceil(len(input_file)/3), 3, sharex=True, sharey=True, figsize=(16, 20))
    axs = axs.ravel()


for counter, target in enumerate(input_file):
    # gets data for plotting
    target_name = target['target_name']
    x_data, y_data = get_data(target['data_file'])
    x_model, y_model = grid_dusty[target['index']]
    y_model = y_model * input_file[counter]['norm']

    # gets supplementary spectra
    if glob("supp_data/*"+target_name.split(' ')[-1]+'*'):
        bonus_spec = []
        for item in os.listdir('./supp_data/'):
            if fnmatch(item, 'IRS*'+target_name.split(' ')[-1]+'*'):
                supp_data = np.array(get_supp_data('supp_data/'+item))
                bonus_spec.append(supp_data)

        # plots supplementary spectra
        for i in range(0, len(bonus_spec)):
            supp_plot_array = Table(np.array([[row[0], row[1]] for row in bonus_spec[i]]), names=('wave', 'lamflam'))
            axs[counter].plot(np.log10(supp_plot_array['wave']), np.log10(
                supp_plot_array['lamflam']), c='k', linewidth=0.5, zorder=2)

        bonus_phot=[]
        for item in os.listdir('./supp_data/'):
            if fnmatch(item, 'phot*'+target_name.split(' ')[-1]+'*'):
                supp_data = np.array(get_supp_data('supp_data/'+item))
                bonus_phot.append(supp_data)

        # plots supplementary phot
        for i in range(0, len(bonus_phot)):
            supp_phot_array = Table(np.array([[row[0], row[1]] for row in bonus_phot[i]]), names=('wave', 'lamflam'))
            axs[counter].scatter(np.log10(supp_phot_array['wave']), np.log10(
                supp_phot_array['lamflam']), facecolor='white', s=20, edgecolor='k', linewidth=0.5, zorder=1)

        bonus_err = []
        for item in os.listdir('./supp_data/'):
            if fnmatch(item, 'err*'+target_name.split(' ')[-1]+'*'):
                supp_data = np.array(get_supp_data('supp_data/'+item))
                bonus_err.append(supp_data)

        # plots supplementary err
        for i in range(0, len(bonus_err)):
            supp_err_array = Table(np.array([[row[0], row[1]] for row in bonus_err[i]]), names=('wave', 'lamflam'))
            supp_err_array['lamflam'][supp_err_array['lamflam'] == 0] = np.nan
            yerror = np.log10(supp_phot_array['lamflam'])*supp_err_array['lamflam']/supp_phot_array['lamflam']
            axs[counter].errorbar(np.log10(supp_phot_array['wave']), np.log10(
                supp_phot_array['lamflam']), yerr=yerror, color='0.3', ls='none', linewidth=0.2, zorder=1)

    # logscale
    x_model = np.log10(x_model)
    y_model = np.log10(y_model)
    x_data = np.log10(x_data)
    y_data = np.log10(y_data)

    # plotting
    if len(input_file) == 1:
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


a = Table.read('fitting_results.csv')



if len(a) == 21:
    additional_data = Table.read('supp_data/GB_additional_data.csv', format='csv')
    a.add_column(additional_data['Av'], index=2)
    a.add_column(additional_data['GB_p'], index=3)
    a.add_column(additional_data['GB_OH'], index=4)

    rgd = Column(((((a['vexp_predicted'])/a['GB_OH']) ** 2) * 200), name='rgd', dtype='int32')
    a.add_column(rgd)

    averages = [('Median'),np.median(a['L']),np.median(a['Av']), np.median(a['GB_p'][a['GB_p'] > 0]), np.median(a['GB_OH']), np.median(a['vexp_predicted']), np.median(a['teff']), np.median(a['tinner']), np.median(a['odep']), np.median(a['mdot']),np.median(a['rgd'])]
    a.add_row(averages)

    a.write('gb_latex_table.csv', format='csv', overwrite=True)
    a['L'] = (a['L']/1000).astype(int)
    a.write('gb_latex_table.txt', format='latex', overwrite=True)

    shutil.copy('gb_latex_table.txt', '/Users/sgoldman/Dropbox/gb_mnras')
    shutil.copy('output_seds.png', '/Users/sgoldman/Dropbox/gb_mnras/figs/GB_seds')