import numpy as np
import matplotlib.pyplot as plt
import pylab
import os
import shutil
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

gb = np.genfromtxt('gb_latex_table.csv', delimiter=',', names=True)
GC = np.genfromtxt('data_GC.csv', delimiter=',', names=True)
LMC = np.genfromtxt('data_LMC.csv', delimiter=',', names=True)

fig = plt.figure(figsize=(4.5, 7))

ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)  # Share y-axes with subplot 1
# Set y-ticks of subplot 2 invisible
plt.setp(ax1.get_xticklabels(), visible=False)


ls = [11128, 30910]
vs = [16, 16]
ms = [2.091E-4, 4.499E-4]

# line
ax1.plot(ls, vs, linestyle='--', c='0.5', linewidth=0.5, zorder=1)

# Plot data
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.yaxis.set_major_formatter(ScalarFormatter())
ax1.yaxis.set_minor_formatter(FormatStrFormatter("%d"))
ax1.yaxis.set_major_formatter(FormatStrFormatter("%d"))
ax1.scatter(gb['L'], gb['GB_OH'], linewidth=0.5, edgecolor='k', c='red', s=16, facecolor='k', marker='o', zorder=4)
ax1.scatter(GC['GC_luminosity'], GC['GC_OH'], linewidth=0.5, edgecolor='k', c='royalblue', s=16,
            facecolor='k', marker='^', zorder=4)
ax1.scatter(LMC['LMC_luminosity'], LMC['LMC_OH'], linewidth=0.5, edgecolor='k', c='limegreen', s=16,
            facecolor='k', marker='s', zorder=4)
ax1.set_xlim(1000, 550000)
ax1.set_ylim(6.01, 31)
ax1.set_ylabel('V$\mathrm{_{exp}}$ (km/s)')
ax1.get_xaxis().set_tick_params(which='both', direction='in')


ax2.set_xscale('log')
ax2.set_yscale('log')
# ax2.xaxis.set_major_formatter(FormatStrFormatter("%d"))
ax2.scatter(gb['L'], gb['mdot'], linewidth=0.5, edgecolor='k', c='red', s=16, facecolor='k',
            label='Galactic Bulge', marker='o')
ax2.scatter(GC['GC_luminosity'], GC['GC_mdot'], linewidth=0.5, edgecolor='k', c='royalblue', s=16,
            facecolor='k', label='Galactic Centre', marker='^')
ax2.scatter(LMC['LMC_luminosity'], LMC['LMC_mdot'], linewidth=0.5, edgecolor='k', c='limegreen', s=16,
            facecolor='k', label='LMC', marker='s')
ax2.set_ylim(0.0000005, 0.006)

# line
ax2.plot(ls, ms, linestyle='--', c='0.5', linewidth=0.5, zorder=1)

ax2.set_xlabel('Luminosity (L${_\\odot}$)')
ax2.set_ylabel('Mass loss (M${_\\odot}$ yr$^{-1}$)')
plt.tight_layout()

leg = plt.legend(loc='lower right')
leg.get_frame().set_edgecolor('0.5')

plt.subplots_adjust(wspace=0.01, hspace=0)


fig = plt.gcf()
fig.set_size_inches(4.5, 7)
fig.savefig('gb_vexp_mdot.png', dpi=300, bbox_inches='tight')
plt.close()

shutil.copy('gb_vexp_mdot.png', '/Users/sgoldman/Dropbox/gb_mnras/figs')
